# study: .sqlite3
# PPM: 20.0
# RT_TOL: 2.0
# MIN_SIGNAL: 100000.0
# MIN_RANGE: 10.0
import sys
import time
import bisect
import sqlite3

start_time = time.time()

batchmode = True

__version__ = "0.7"


#
# Based on: https://github.com/suzaku/plain_obj
#


def make_constructor(fields):
    assignments = "\n".join([f"    self.{f} = {f}" for f in fields])
    parameter_lists = ", ".join(fields)
    source = "def __init__(self, %s):\n%s" % (parameter_lists, assignments)
    namespace = {}
    exec(source, namespace)
    return namespace["__init__"]


def make_lt(key_field):
    source = (
        f"def __lt__(self, other):\n    return self.{key_field} < other.{key_field}"
    )
    namespace = {}
    exec(source, namespace)
    return namespace["__lt__"]


def sortable(type_name, field_names):
    if isinstance(field_names, str):
        # names separated by whitespace and/or commas
        field_names = field_names.replace(",", " ").split()
    return type(
        type_name,
        (object,),
        {
            "__slots__": field_names,
            "__init__": make_constructor(field_names),
            "__lt__": make_lt(field_names[0]),
        },
    )


ppm = PPM / 1000000.0


MS1_DBNAME = study

#
# Set up our cursor
#
con = sqlite3.connect(MS1_DBNAME)
cur = con.cursor()

#
# First, a sanity check (somewhat repetitive to later code, but this is temporary and will be removed when all unsafe raw2sql
# generation will have been eliminated from the lab, along with any incorrect sqlite3 file...
#

sample_num = cur.execute(
    "SELECT COUNT (DISTINCT name) from rawfile where ID > 0"
).fetchone()[0]
observed = cur.execute(
    "SELECT COUNT (DISTINCT rawfile) from scans where rawfile > 0"
).fetchone()[0]

if sample_num != observed:
    print("Sanity check failure: expected samples != observed files!!!", file=sys.stderr, flush=True)
    sys.exit(-1)

#
# Second, load mass_translation and time_translation factors, if they are available...
# TODO: Potentially factor this code out into a module that manages access to our db format. 
#

try:
    mass_translation_factor = cur.execute(
        "SELECT value from sequence where attribute = 'mass_translation_factor'"
    ).fetchone()[0]
    time_translation_factor = cur.execute(
        "SELECT value from sequence where attribute = 'time_translation_factor'"
    ).fetchone()[0]
except:
    mass_translation_factor = 10000  # 1 = 0.0001 Da
    time_translation_factor = 1000   # 1 = 0.001 seconds


PEAK_SQL = """SELECT ms1_peaks.intensity, ms1_peaks.mz, ms1_peaks.rt, ms1_peaks.rawfile FROM ms1_peaks WHERE ms1_peaks.rawfile > 0 AND ms1_peaks.intensity > ? ORDER BY ms1_peaks.intensity DESC"""

print(f"Starting the mega-query...", file=sys.stderr, flush=True)
start_query = time.time()
all_peaks = list(cur.execute(PEAK_SQL, (MIN_SIGNAL,)))
stop_query = time.time()

print(f"Finished mega-query in {stop_query-start_query:.2f} seconds.", file=sys.stderr, flush=True)

print(f"The total number of peaks is: {len(all_peaks)}", file=sys.stderr, flush=True)

print(f"Starting the peak processing...", file=sys.stderr, flush=True)
start_process = time.time()

Slot = sortable("Slot", "mz rt file min max")
slots = []
rt_max = 0
for (intensity, mz, rt, rfile) in all_peaks:
    # NOTE: rt and mz are still raw integers unmodified by time and mass factors!!!
    mz = float(mz / mass_translation_factor)
    rt = float(rt / (60 * time_translation_factor))        
    potential_slot = Slot(mz, rt, rfile, intensity, intensity)
    if rt > rt_max:
        rt_max = rt
    if mz > 0.0:
        lower = (1.0 - ppm) * mz
        upper = (1.0 + ppm) * mz
    else:
        lower = (1.0 + ppm) * mz
        upper = (1.0 - ppm) * mz
    potential_slot.mz = lower
    lower_index = bisect.bisect_left(slots, potential_slot)
    potential_slot.mz = upper
    upper_index = bisect.bisect_right(slots, potential_slot)
    potential_slot.mz = mz
    assigned = False
    if (lower_index != len(slots)) and (upper_index > 0):
        for offset in range(lower_index, upper_index):
            slot = slots[offset]
            if abs(slot.rt - rt) < RT_TOL:
                assigned = True
                if slot.file == rfile and intensity < slot.min:
                    slot.min = intensity
    if not assigned:
        bisect.insort(slots, potential_slot)

stop_process = time.time()
print(
    f"Finished peak processing {len(all_peaks)} peaks into {len(slots)} potential features in {stop_process-start_process:.2f} seconds.",
    file=sys.stderr, flush=True
)

out = open(__file__[:-3] + ".features", "w")
print(
    "Metabolite\tFormula\tIon Type\tRT (min)\tPercentile",
    file=out,
)

total_feature_counter = 0
for slot in slots:
    if slot.max / slot.min > MIN_RANGE:
        total_feature_counter += 1

final_feature_counter = total_feature_counter
for slot in slots:
    if slot.max / slot.min > MIN_RANGE:
        final_feature_counter -= 1
        polarity = "+"
        mz_name = abs(slot.mz)
        if slot.mz < 0:
            polarity = "-"
        print(
            f"Feature_{polarity}_{mz_name:.4f}_{slot.rt:.1f}\t{slot.mz :.4f}\t\t{slot.rt:.3f}\t{final_feature_counter / total_feature_counter:0.3f}",
            file=out,
        )
out.close()

stop_time = time.time()

print(
    f"Ungrid processed {len(all_peaks)} peaks into {total_feature_counter} features in {stop_time - start_time :.2f} seconds.",
    file=sys.stderr, flush=True
)
