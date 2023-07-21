# study: .sqlite3
# neutral_loss: 132.0420
# abs_mz: False True
# PPM: 20.0
# RT_TOL: 0.5
# output: .features
import sys
import time
import bisect
import sqlite3

start_time = time.time()

__version__ = "0.5"


class Thing:
    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)


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


NLOSS = neutral_loss

if NLOSS > 0:
    mz_low = NLOSS * (1.0 - ppm)
    mz_high = NLOSS * (1.0 + ppm)
else:
    mz_low = NLOSS * (1.0 + ppm)
    mz_high = NLOSS * (1.0 - ppm)
mz_low = round(mz_low * mass_translation_factor)
mz_high = round(mz_high * mass_translation_factor)

if abs_mz == "True":
    PEAK_SQL = f"SELECT ms2_peaks.precursor, ms2_peaks.mz, ms2_peaks.rt, ms2_peaks.rawfile FROM ms2_peaks WHERE abs(ms2_peaks.precursor - ms2_peaks.mz) BETWEEN ? AND ? ORDER BY intensity DESC"

else:
    PEAK_SQL = f"SELECT ms2_peaks.precursor, ms2_peaks.mz, ms2_peaks.rt, ms2_peaks.rawfile FROM ms2_peaks WHERE ms2_peaks.precursor - ms2_peaks.mz BETWEEN ? AND ? ORDER BY intensity DESC"

print(f"Starting the neutral loss search...", file=sys.stderr, flush=True)
start_query = time.time()
all_peaks = list(cur.execute(PEAK_SQL, (mz_low, mz_high)))
stop_query = time.time()

print(
    f"Finished precursor-query in {stop_query-start_query:.2f} seconds.",
    file=sys.stderr, flush=True
)

print(
    f"The total number of peaks with neutral loss is: {len(all_peaks)}", file=sys.stderr, flush=True
)

print(f"Outputing neutral losses", file=sys.stderr, flush=True)
start_process = time.time()

Slot = sortable("Slot", "precursor mz rt file")
slots = []
rt_max = 0
for (precursor, mz, rt, rfile) in all_peaks:
    precursor = float(precursor / mass_translation_factor)
    mz = float(mz / mass_translation_factor)
    rt = float(rt / (60 * time_translation_factor))
    potential_slot = Slot(precursor, mz, rt, rfile)
    if rt > rt_max:
        rt_max = rt
    if precursor > 0.0:
        lower = (1.0 - ppm) * precursor
        upper = (1.0 + ppm) * precursor
    else:
        lower = (1.0 + ppm) * precursor
        upper = (1.0 - ppm) * precursor
    potential_slot.precursor = lower
    lower_index = bisect.bisect_left(slots, potential_slot)
    potential_slot.precursor = upper
    upper_index = bisect.bisect_right(slots, potential_slot)
    potential_slot.precursor = precursor
    assigned = False
    if (lower_index != len(slots)) and (upper_index > 0):
        for offset in range(lower_index, upper_index):
            slot = slots[offset]
            if abs(slot.rt - rt) < RT_TOL:
                assigned = True
    if not assigned:
        bisect.insort(slots, potential_slot)

stop_process = time.time()
print(
    f"Finished processing {len(all_peaks)} neutral losses into {len(slots)} potential parents in {stop_process-start_process:.2f} seconds.",
    file=sys.stderr, flush=True
)
out = open(__file__[:-3] + ".features", "w")
print(
    "Metabolite\tFormula\tIon Type\tRT Start (min)\tRT End (min)\tm/z Tolerance (ppm)\tRT Tolerance (min)",
    file=out,
)

for slot in slots:
    rt_start = max(0, slot.rt - RT_TOL)
    rt_stop = min(slot.rt + RT_TOL, rt_max)
    print(
        f"Parent_{slot.precursor :.4f}_{NLOSS :.4f}_{slot.mz :.4f}_{slot.rt :.1f}\t{slot.precursor:.4f}\t\t{rt_start :.2f}\t{rt_stop :.2f}\t{PPM}\t{RT_TOL}",
        file=out,
    )
out.close()
