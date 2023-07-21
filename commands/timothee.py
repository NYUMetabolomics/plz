# study: .sqlite3
import sys
import os
import time
import subprocess
import sqlite3
import lib2nist

start_time = time.time()

__version__ = "0.1.1"

# class Thing:
#     def __getitem__(self, key):
#         return getattr(self, key)

#     def __setitem__(self, key, value):
#         setattr(self, key, value)


#
# Set up our cursor
#
con = sqlite3.connect(study)
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



FILE_SQL = """SELECT rawfile.id, rawfile.name FROM rawfile WHERE rawfile.id > 0 ORDER BY rawfile.id ASC"""
SCAN_SQL = """SELECT scans.scan_ID, scans.rt, scans.precursor, scans.scan_type FROM scans WHERE scans.rawfile = ? ORDER BY scans.scan_ID ASC"""
PEAK_SQL = """SELECT ms2_peaks.rawfile, ms2_peaks.rt, ms2_peaks.mz, ms2_peaks.intensity FROM ms2_peaks WHERE ms2_peaks.rawfile > 0 ORDER BY ms2_peaks.rawfile ASC, ms2_peaks.rt ASC, ms2_peaks.mz ASC"""

all_filenames = {}
for x in cur.execute(FILE_SQL):
    all_filenames[x[0]] = x[1]

seq_name = study[:-8]

scan_loader_start = time.time()
scan_loader_counter = 0
scan_loader_scan_counter = 0
all_scans = {}
for (rawid, fname) in all_filenames.items():
    scan_loader_counter += 1
    print(f"processing {fname} ({scan_loader_counter}) / ({len(all_filenames.items())})", file=sys.stderr, flush=True)
    scans = {}
    prev_ms1 = None
    for x in cur.execute(SCAN_SQL, (rawid,)):
        if x[3] == "MS1":
            prev_ms1 = x[1]
        else:
            scan_loader_scan_counter += 1
            scans[x[1]] = [x[0], x[2], prev_ms1]
    all_scans[rawid] = scans
scan_loader_stop = time.time()
print(f"{scan_loader_scan_counter} scans loaded in {scan_loader_stop - scan_loader_start:.2f} seconds", file=sys.stderr, flush=True)

#
# MGF Generation and Search
#

mgf_make_start = time.time()

mgf = open(seq_name + ".mgf", 'w')

start_query = time.time()
prev_rt = -1
prev_rawid = None
scan_count = 0
peak_count = 0
peaks = []  # Technically this gets taken care of in the loop, but still, defensive programming etc...
peak_counter = 0  # TODO: duplicate counter used purely for visual feedback, should be removed.
for (rawid, rt, mz, intensity) in cur.execute(PEAK_SQL):
    # NOTE: rt and mz are still raw integers unmodified by time and mass factors!!!
    peak_counter += 1
    if not (peak_counter % 100000):
        print(f"peaks processed = {peak_counter} ", flush=True)
    if rawid != prev_rawid or rt != prev_rt:
        # Finalize Scan...
        if scan_count > 0:
            for (m, i) in peaks:
                print(f"{m} {i}", file=mgf)
            print("END IONS\n", file=mgf)
        peaks = []
        prev_rt = rt
        prev_rawid = rawid
        scan_count += 1
        pepmass = all_scans[rawid][rt][1]
        print("BEGIN IONS", file=mgf)
        if pepmass > 0:
            print(f"TITLE={all_filenames[rawid]}.{all_scans[rawid][rt][0]}.+", file=mgf)
            print("CHARGE=127+", file=mgf)
        else:
            print(f"TITLE={all_filenames[rawid]}.{all_scans[rawid][rt][0]}.-", file=mgf)
            print("CHARGE=128-", file=mgf)
        print(f"RTINSECONDS={float(rt / time_translation_factor)}", file=mgf)
        if pepmass > 0:
            print(f"PEPMASS={float(pepmass / mass_translation_factor)}", file=mgf)
        else:
            print(f"PEPMASS={-1 * float(pepmass / mass_translation_factor)}", file=mgf)
    peak_count += 1
    if mz < 0:
        peaks = [(float( (-mz) / mass_translation_factor), intensity)] + peaks
    else:
        peaks.append((float(mz / mass_translation_factor), intensity))
# Finalize Last Scan...
if scan_count > 0:
    for (m, i) in peaks:
        print(f"{m} {i}", file=mgf)
    print("END IONS\n", file=mgf)
stop_query = time.time()
mgf.close()

mgf_make_stop = time.time()
print(f"MGF created in {mgf_make_stop - mgf_make_start:.2f} seconds.", file=sys.stderr, flush=True)
print(f"sql2mgf phase processed {peak_count} peaks and {scan_count} scans from {len(all_filenames)} files in {mgf_make_stop - start_time :.2f} seconds.", file=sys.stderr, flush=True)
print("", file=sys.stderr, flush=True)
print("----------------", file=sys.stderr, flush=True)


msp_name = f"{seq_name + '.MSP'}"
msp = open(msp_name, 'w')
name = None
precursor = None
current_id = 1
peaks = []
with open(seq_name + ".mgf") as mgf:
    for line in mgf:
        if line.startswith("END IONS"):
            print(f"Name: {name}", file=msp)
            print(f"PrecursorMZ: {precursor}", file=msp)
            print(f"ID: {current_id}", file=msp)
            print(f"Num peaks: {len(peaks)}", file=msp)
            for peak in peaks:
                print(peak, file=msp)
            print("", file=msp)
            name = None
            precursor = None
            current_id += 1
            peaks = []
            if not current_id % 65536:
                # This allows us to use the default KEEP_IDS setting in LIB2NIST and nevertheless
                # circumvent its known bug regarding IDS at multiples of 65536.
                current_id += 1
            continue
        if line.startswith("TITLE="):
            name = line.strip()[7:]
            continue
        if line.startswith("RTINSECONDS="):
            name += "." + (line.strip()[12:].replace(".",""))
            continue
        if line.startswith("PEPMASS="):
            precursor = line.strip()[8:]
            continue
        if line[0] in "0123456789":
            peaks.append(line.strip())
            continue
msp.close()

if lib2nist.platform == "win32":
    # print([lib2nist.path, msp_name, f"{os.getcwd()}\\", "/AccuratePeakMZ", "/MsmsOnly:Y"])
    # print(f"{lib2nist.path} {msp_name} {os.getcwd()}\\ /AccuratePeakMZ /MsmsOnly:Y")
    subprocess.run([lib2nist.path, msp_name, f"{os.getcwd()}\\", "/AccuratePeakMZ", "/MsmsOnly:Y"], shell=False)
    # os.system(f"{lib2nist.path} {msp_name} {os.getcwd()}\\ /AccuratePeakMZ /MsmsOnly:Y")
else:
    print("timothee currently unavailable on non-win32 systems!!!", file=sys.stderr, flush=True)
    sys.exit(-1)
