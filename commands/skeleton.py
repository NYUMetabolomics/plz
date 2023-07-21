# study: .sqlite3
# features: .features
# MZ_TOLERANCE: 15.0
# RT_TOLERANCE: 0.2
# RT_WINDOW: 0.5
# output: .quantified
import sys
import time
import sqlite3
import actions
import formula_actions
import nist_ion_descriptions
import formula
import columns

# Used to produce these here:
#
# \tRT Start (min)\tRT End (min)\tm/z Tolerance (ppm)\tRT Tolerance (min)
# \t{max(rt - RT_WINDOW, 0.0):.02f}\t{rt + RT_WINDOW:.02f}\t{MZ_TOLERANCE:.2f}\t{RT_TOLERANCE:.2f}
#
# Using these:
#
# MZ_TOLERANCE: 15.0
# RT_TOLERANCE: 0.2
# RT_WINDOW: 0.25
#
# MZ_TOLERANCE is specified in ppm and corresponds to the tolerance parameter ultimately used by skeleton quantitation.
# RT_TOLERANCE is specified in minutes and corresponds to the tolerance used by the skeleton quantitation.
# RT_WINDOW is used to bracket the location of the MS2 ID -- it indirectly controls first phase of skeleton.

start_time = time.time()

__version__ = "1.1.0"

MS1_DBNAME = study
input_filename = features

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


#
# Load immutable scan information, which will be necessary (e.g. for default zero intensity per rawfile x rt)
#

SAMPLES_SQL = """SELECT DISTINCT name from rawfile where ID > 0 order by ID"""

samples = [entry[0] for entry in cur.execute(SAMPLES_SQL)]
# NOTE: sample[0] often corresponds to ID = 1 -- IDs start at 1 so offsets in samples will _not_ immediately correspond to ID!!!
SCAN_SQL = """SELECT scans.rawfile, scans.rt from scans where scans.rawfile > 0 and scans.polarity = ? and scans.scan_type = 'MS1' order by scans.rawfile, scans.rt"""

pos_rawfile_rt_pairs = list(cur.execute(SCAN_SQL, ("+")))
neg_rawfile_rt_pairs = list(cur.execute(SCAN_SQL, ("-")))

EXICS_SQL = """SELECT ms1_peaks.rawfile, ms1_peaks.rt, ms1_peaks.intensity, ms1_peaks.mz from ms1_peaks where ms1_peaks.rawfile > 0 and ms1_peaks.mz >= ? and ms1_peaks.mz <= ? and ms1_peaks.rt >= ? and ms1_peaks.rt <= ? order by ms1_peaks.rawfile, ms1_peaks.rt, ms1_peaks.intensity"""

# cur.close()


def exics(polarity, mz, rt_start, rt_stop, mz_tol):
    # retention times are assumed to be provided in minutes...
    # TODO: modify signature to remove polarity, for now keep in mind that mz is
    #       already negative when polarity is "-"
    # cur = con.cursor()
    xics = {}
    if polarity == "+":
        rawfile_rts = pos_rawfile_rt_pairs
    else:
        rawfile_rts = neg_rawfile_rt_pairs
    for (rawfile, rt) in rawfile_rts:
        if rawfile not in xics:
            xics[rawfile] = {}
        # Note that while the xics include data-pairs for all RTs,
        # the ones with non-zero data are, in fact, constrained by
        # rt_start and rt_stop...
        xics[rawfile][rt] = (0, 0)  # Setting the intensity to zero could mask a "real" zero in the data (for negative m/z values)... Assuming those are allowed.

    if polarity == "+":
        mz_low = mz * (1.0 - mz_tol)
        mz_high = mz * (1.0 + mz_tol)
    else:
        mz_high = mz * (1.0 - mz_tol)
        mz_low = mz * (1.0 + mz_tol)

    mz_low = round(mz_low * mass_translation_factor)
    mz_high = round(mz_high * mass_translation_factor)
    rt_start = round(rt_start * 60 * time_translation_factor)  # time_translation_factor is in seconds not minutes
    rt_stop = round(rt_stop * 60 * time_translation_factor)    # time_translation_factor is in seconds not minutes

    cur.execute(EXICS_SQL, (mz_low, mz_high, rt_start, rt_stop))
    for (rawfile, rt, intensity, omz) in cur:
        if polarity == "+":
            xics[rawfile][rt] = max(xics[rawfile][rt], (intensity, omz))
        else:
            xics[rawfile][rt] = max(xics[rawfile][rt], (intensity, -1 * omz))  # Here a (0, -100) would lose out to the default (0, 0)

    float_min_time_factor = 60.0 * float(time_translation_factor)
    float_mass_factor = float(mass_translation_factor)
    xic_list = []
    prev_rawfile = 0
    for rawfile in xics:
        assert (
            rawfile > prev_rawfile
        )  # this should be true based on dictionary semantics in python, but better safe than sorry :-)
        prev_rawfile = rawfile
        xic_list.append(
            list(
                map(
                    lambda rt__int_mz: (
                        rt__int_mz[0] / float_min_time_factor,
                        (rt__int_mz[1][0], rt__int_mz[1][1] / float_mass_factor),
                    ),
                    xics[rawfile].items(),
                )
            )
        )

    # cur.close()
    return xic_list


# reference_masses = [{"name": "Hydrogen", "symbol": "H", "mass": 1007825035}, {"name": "Silicon", "symbol": "Si", "mass": 27976926530}, {"name": "Lithium", "symbol": "Li", "mass": 7016003000}, {"name": "Boron", "symbol": "B", "mass": 11009305500}, {"name": "Carbon", "symbol": "C", "mass": 12000000000}, {"name": "Nitrogen", "symbol": "N", "mass": 14003074000}, {"name": "Oxygen", "symbol": "O", "mass": 15994914630}, {"name": "Fluorine", "symbol": "F", "mass": 18998403220}, {"name": "Sodium", "symbol": "Na", "mass": 22989767700}, {"name": "Magnesium", "symbol": "Mg", "mass": 23985042300}, {"name": "Phosphorous", "symbol": "P", "mass": 30973762000}, {"name": "Sulfur", "symbol": "S", "mass": 31972070700}, {"name": "Chlorine", "symbol": "Cl", "mass": 34968852720}, {"name": "Potassium", "symbol": "K", "mass": 38963707400}, {"name": "Calcium", "symbol": "Ca", "mass": 39962590600}, {"name": "Chromium", "symbol": "Cr", "mass": 51940509800}, {"name": "Manganese", "symbol": "Mn", "mass": 54938047100}, {"name": "Iron", "symbol": "Fe", "mass": 55934939300}, {"name": "Nickel", "symbol": "Ni", "mass": 57935346200}, {"name": "Cobalt", "symbol": "Co", "mass": 58933197600}, {"name": "Copper", "symbol": "Cu", "mass": 62929598900}, {"name": "Zinc", "symbol": "Zn", "mass": 63929144800}, {"name": "Arsenic", "symbol": "As", "mass": 74921594200}, {"name": "Bromine", "symbol": "Br", "mass": 78918336100}, {"name": "Selenium", "symbol": "Se", "mass": 79916519600}, {"name": "Molybdenum", "symbol": "Mo", "mass": 97905407300}, {"name": "Palladium", "symbol": "Pd", "mass": 105903478000}, {"name": "Silver", "symbol": "Ag", "mass": 106905092000}, {"name": "Cadmium", "symbol": "Cd", "mass": 113903357000}, {"name": "Iodine", "symbol": "I", "mass": 126904473000}, {"name": "Gold", "symbol": "Au", "mass": 196966543000}, {"name": "Mercury", "symbol": "Hg", "mass": 201970617000}]

# atomic_mass = {}
# symbols = []
# for entry in reference_masses:
#     atomic_mass[entry["symbol"]] = entry["mass"] * 0.000000001
#     symbols.append(entry["symbol"])

# symbols.sort(key=len, reverse=True)

proton = 1.00727647
electron = 0.00054858


C13_delta = 1.003354838  # https://en.wikipedia.org/wiki/Isotopes_of_carbon
N15_delta = 0.9970348934  # https://en.wikipedia.org/wiki/Isotopes_of_nitrogen
O18_delta = 2.0042463804  # https://en.wikipedia.org/wiki/Isotopes_of_oxygen

digits = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]


def find_first_nonmember(string, reference):
    i = 0
    while i < len(string):
        if string[i] not in reference:
            break
        i += 1
    return i


# def positive_integer_prefix(string):
#     return find_first_nonmember(string, digits)


# def parse_formula(formula):
#     complete = formula[:]
#     mass = 0.0
#     c_counter = 0
#     n_counter = 0
#     o_counter = 0
#     while formula != "":
#         key = ""
#         for symbol in symbols:
#             if formula.startswith(symbol):
#                 key = symbol
#                 break
#         if key == "":
#             return ["The formula fragment: '" +
#                     formula + "' of the formula '" +
#                     complete + "' is missing a symbol prefix...", -1]

#         formula = formula[len(key):]
#         prefix = positive_integer_prefix(formula)
#         if prefix == 0:
#             count = 1
#         else:
#             count = int(formula[0:prefix])
#             formula = formula[prefix:]

#         submass = count * atomic_mass[key]
#         if key == "C":
#             c_counter += count
#         if key == "N":
#             n_counter += count
#         if key == "O":
#             o_counter += count
#         mass += submass
#     return (mass, c_counter, n_counter, o_counter)


def parse_labeling_part(
    metabolite,
    labeling_part,
    is_molecular_formula,
    c_counter,
    n_counter,
    o_counter,
    part_of_pair=False,
):
    complete_labeling_part = labeling_part[:]
    M_label = ""
    if (
        not labeling_part.startswith("13C")
        and not labeling_part.startswith("15N")
        and not labeling_part.startswith("18O")
    ):
        print(
            f"Invalid labeling label {labeling_part} for metabolite {metabolite}!!!",
            file=sys.stderr, flush=True
        )
        sys.exit(-1)
    if labeling_part.startswith("13C"):
        M_label = "13C"
    if labeling_part.startswith("15N"):
        M_label = "15N"
    if labeling_part.startswith("18O"):
        M_label = "18O"
    labeling_part = labeling_part[
        3:
    ]  # not sure what this means, may need to change for 18O
    if (not labeling_part.startswith("[")) and (not is_molecular_formula):
        print(
            f"Cannot expand {complete_labeling_part} for {metabolite} when the provided 'formula' is an m/z value!!!",
            file=sys.stderr, flush=True
        )
        sys.exit(-1)
    if labeling_part.startswith("["):
        labeling_part = labeling_part[1:(-1)]
        if "-" in labeling_part:
            range_spec = labeling_part.split("-")
            M_range = list(range(int(range_spec[0]), int(range_spec[1]) + 1))
        else:
            M_range = [int(labeling_part)]
    else:
        if M_label == "13C":
            if part_of_pair:
                M_range = [c_counter]
            else:
                M_range = list(range(0, c_counter + 1))
        else:
            if M_label == "15N":
                if part_of_pair:
                    M_range = [n_counter]
                else:
                    M_range = list(range(0, n_counter + 1))
            else:
                if M_label == "18O":
                    if part_of_pair:
                        M_range = [o_counter]
                    else:
                        M_range = list(range(0, o_counter + 1))
                else:
                    print("Unknown Label", file=sys.stderr, flush=True)
    return (M_label, M_range)


def calc_labeling_label_delta(M_label, M_value):
    if M_label == "13C":
        return M_value * C13_delta
    if M_label == "15N":
        return M_value * N15_delta
    if M_label == "18O":
        return M_value * O18_delta
    print(
        "Incorrect labeling requirements encountered: {M_label}, {M_value}",
        file=sys.stderr, flush=True
    )
    sys.exit(-1)


def parse_labeling_requirements(
    metabolite, labeling, base_mz, is_molecular_formula, c_counter, n_counter, o_counter
):
    M_suffixes = []
    M_mzs = []

    if labeling.startswith('"'):
        labeling = labeling[1:(-1)]

    if not labeling:  # An empty labeling field will be ignored
        return ([], [])

    if "," in labeling:
        labeling_parts = labeling.split(",")
        (M1_label, M1_range) = parse_labeling_part(
            metabolite,
            labeling_parts[0],
            is_molecular_formula,
            c_counter,
            n_counter,
            o_counter,
            True,
        )
        (M2_label, M2_range) = parse_labeling_part(
            metabolite,
            labeling_parts[1],
            is_molecular_formula,
            c_counter,
            n_counter,
            o_counter,
            True,
        )
        assert len(M1_range) == 1
        assert len(M2_range) == 1
        M_suffixes = [f"-{M1_label},{M2_label}-{M1_range[0]},{M2_range[0]}"]
        if base_mz > 0:
            M_mzs = [
                base_mz
                + calc_labeling_label_delta(M1_label, M1_range[0])
                + calc_labeling_label_delta(M2_label, M2_range[0])
            ]
        else:
            M_mzs = [
                base_mz
                - calc_labeling_label_delta(M1_label, M1_range[0])
                - calc_labeling_label_delta(M2_label, M2_range[0])
            ]
    else:
        (M_label, M_range) = parse_labeling_part(
            metabolite, labeling, is_molecular_formula, c_counter, n_counter, o_counter
        )
        for M_value in M_range:
            M_suffixes.append(f"-{M_label}-{M_value}")
            if base_mz > 0:
                M_mzs.append(base_mz + calc_labeling_label_delta(M_label, M_value))
            else:
                M_mzs.append(base_mz - calc_labeling_label_delta(M_label, M_value))
    return (M_suffixes, M_mzs)


metabolites = set()

metabolite_2_formula = {}
metabolite_2_mz = {}
metabolite_2_observed = {}
metabolite_2_rs = {}
metabolite_2_measurements = {}

print(f"loading XIC request file: {input_filename}", file=sys.stderr, flush=True)

(rows, unmatched, unexpected) = columns.loader(
    input_filename,
    {
        "Metabolite": {
            "field": "metabolite",
            "constructor": str,
            "required": True,
            "description": "Metabolite name.",
        },
        "Formula": {
            "field": "formula",
            "constructor": str,
            "required": True,
            "description": "Formula or m/z value.",
        },
        "InChIKey": {
            "field": "inchikey",
            "constructor": str,
            "required": False,
            "description": "InChIKey.",
        },
        #  "Polarity (z)": {
        #     "field": "polarity",
        #     "constructor": enum("+ -"),
        #     "required": True,
        #     "description": "Polarity: (+/-)."
        # },
        "Labeling": {
            "field": "labeling",
            "constructor": str,
            "required": False,
            "description": "Labeling scheme assumed for the given metabolite.",
        },
        "Ion Type": {
            "field": "ion_type",
            "constructor": str,
            "required": True,
            "description": "Ionization type.",
        },
        # Either RT or RT Start+End must be present...
        "RT (min)": {
            "field": "rt",
            "constructor": float,
            "required": False,
            "description": "RT (in minutes).",
        },
        "RT Start (min)": {
            "field": "rt_start",
            "constructor": float,
            "required": False,
            "description": "RT window start (in minutes).",
        },
        "RT End (min)": {
            "field": "rt_stop",
            "constructor": float,
            "required": False,
            "description": "RT window end (in minutes).",
        },
        "m/z Tolerance (ppm)": {
            "field": "mz_tol",
            "constructor": float,
            "required": False,
            "description": "m/z tolerance in PPM.",
        },
        "RT Tolerance (min)": {
            "field": "rt_tol",
            "constructor": float,
            "required": False,
            "description": "RT tolerance (in minutes).",
        },
        "FDR": {
            "field": "fdr",
            "constructor": int,
            "required": False,
            "description": "FDR of Metabolite ID.",
        },
    },
)


def refined_max(an_xic, tight_start, tight_stop):
    return max(
        (intensity, mz, t)
        for (t, (intensity, mz)) in an_xic
        if tight_start <= t <= tight_stop
    )  # This will fail if tight_start to tight_stop is too narrow to afford even one scan per sample!


def imposed_max(an_xic, rt):
    for (t, (intensity, mz)) in an_xic:
        if t == rt:
            return (intensity, mz, t)
    print("imposed rt did not exist in the target xic!", file=sys.stderr, flush=True)
    sys.exit(-1)


def robust_max(aList):
    if not aList:
        return 0.0
    return max(aList)


output = open(__file__[:-3] + ".quantified", "w")
if "inchikey" in unmatched and "labeling" in unmatched and "fdr" in unmatched:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )
elif "inchikey" in unmatched and "labeling" in unmatched and "fdr" not in unmatched:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "FDR",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )
elif "inchikey" in unmatched and "fdr" in unmatched:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "Labeling",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "is_global_winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )
elif "inchikey" in unmatched and "fdr" not in unmatched:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "FDR",
                "Labeling",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "is_global_winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )
elif "labeling" in unmatched and "fdr" in unmatched:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "InChIKey",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )
elif "labeling" in unmatched and "fdr" not in unmatched:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "InChIKey",
                "FDR",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )
elif "fdr" in unmatched:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "InChIKey",
                "Labeling",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "is_global_winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )
else:
    print(
        "\t".join(
            [
                "Metabolite",
                "Formula",
                "InChIKey",
                "FDR",
                "Labeling",
                "Ion Type",
                "RT Start (min)",
                "RT End (min)",
                "m/z Tolerance (ppm)",
                "RT Tolerance (min)",
                "mz",
                "obs_mz",
                "ppm",
                "winner",
                "is_global_winner",
                "RT",
                "RT_min",
                "RT_max",
                "RT_range",
                "detections",
            ]
            + samples
        ),
        file=output,
    )


def process_entry(an_entry):
    (row, M_suffixes, M_mzs, M_measurements) = an_entry
    mega_max = []
    for M_offset in range(len(M_suffixes)):
        measurements = M_measurements[M_offset]
        maxima = list(
            map(lambda x: robust_max([pair for (rt, pair) in x]), measurements)
        )
        winner_pair = max(enumerate(maxima), key=lambda i_max_int: i_max_int[1])  # winner_pair = (offset, (intensity, mz))
        # winner_sample = samples[winner_pair[0]]
        (the_rt, the_pair) = max(measurements[winner_pair[0]], key=lambda x: x[1])  # the_pair should equal winner_pair[1]
        mega_max.append(
            (the_pair[0], M_offset, the_rt)  # the_pair[0] == intensity
        )  # (what intensity, which isotope, when)
    mega_max.sort(reverse=True)
    mega_winner_intensity = mega_max[0][0]
    # mega_winner_sample = mega_max[0][1]
    mega_winner_offset = mega_max[0][1]
    mega_winner_rt = mega_max[0][2]
    # mega_winner_mz = mega_max[0][3]
    sample_winner_rts = []
    for sample_offset in range(len(M_measurements[0])):
        isotope_maxima = []
        for M_offset in range(len(M_suffixes)):
            try:
                isotope_maxima.append(
                    refined_max(
                        M_measurements[M_offset][sample_offset],
                        mega_winner_rt - row.rt_tol,
                        mega_winner_rt + row.rt_tol,
                    )
                )
            except:
                print(f"--> skipping row with insufficient coverage in sample {samples[sample_offset]} (in range: [{mega_winner_rt - row.rt_tol} - {mega_winner_rt + row.rt_tol}])!!!", file=sys.stderr, flush=True)
                return
        sample_winner = max(isotope_maxima)
        sample_winner_rts.append(sample_winner[2])
    for M_offset in range(len(M_suffixes)):
        measurements = M_measurements[M_offset]
        if M_offset == mega_winner_offset:
            is_global_winner = "Yes"
        else:
            is_global_winner = "No"
        if mega_winner_intensity == 0.0:
            rt_min = 0.0
            rt_max = 0.0
            rt_range = 0.0
            finalized = ["0"] * len(measurements)
            local_winner_mz = 0.0
            local_winner_sample = ""
            local_winner_rt = 0.0
            local_max = 0.0
            ppm = 0.0
            detections = 0
        else:
            # finalized_vals = list(map(lambda x: refined_max(x, mega_winner_rt - row.rt_tol, mega_winner_rt + row.rt_tol), measurements))
            finalized_vals = list(
                map(lambda x: imposed_max(*x), zip(measurements, sample_winner_rts))
            )
            local_winner_mz = 0.0
            local_winner_sample = ""
            local_winner_rt = 0.0
            local_max = 0.0
            sample_offset = -1
            rts = []
            for x in finalized_vals:
                sample_offset += 1
                if x[0] > 0.0:
                    rts.append(x[2])
                if x[0] > local_max:
                    local_max = x[0]
                    local_winner_mz = x[1]
                    local_winner_rt = x[2]
                    local_winner_sample = samples[sample_offset]
            detections = len(rts)
            if not rts:
                rt_min = 0
                rt_max = 0
                rt_range = 0
                ppm = 0
            else:
                rt_min = min(rts)
                rt_max = max(rts)
                rt_range = rt_max - rt_min
                if row.polarity == "+":
                    ppm = (
                        1000000.0
                        * (local_winner_mz - M_mzs[M_offset])
                        / M_mzs[M_offset]
                    )
                else:
                    ppm = (
                        1000000.0
                        * (((-1) * local_winner_mz) - M_mzs[M_offset])
                        / M_mzs[M_offset]
                    )

            finalized = list(
                map(lambda x: "%.0f" % x, map(lambda pair: pair[0], finalized_vals))
            )
        if "inchikey" in unmatched and "labeling" in unmatched and "fdr" in unmatched:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formula + M_suffixes[M_offset],
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )
        elif "inchikey" in unmatched and "labeling" in unmatched and "fdr" not in unmatched:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formula + M_suffixes[M_offset],
                        "%d" % row.fdr,
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )            
        elif "inchikey" in unmatched and "fdr" in unmatched:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formula + M_suffixes[M_offset],
                        row.labeling,
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        is_global_winner,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )
        elif "inchikey" in unmatched and "fdr" not in unmatched:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formula + M_suffixes[M_offset],
                        "%d" % row.fdr,
                        row.labeling,
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        is_global_winner,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )
        elif "labeling" in unmatched and "fdr" in unmatched:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formula + M_suffixes[M_offset],
                        row.inchikey + M_suffixes[M_offset],
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )
        elif "labeling" in unmatched and "fdr" not in unmatched:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formula + M_suffixes[M_offset],
                        row.inchikey + M_suffixes[M_offset],
                        "%d" % row.fdr,
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )
        elif "fdr" in unmatched:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formulas + M_suffixes[M_offset],
                        row.inchikey + M_suffixes[M_offset],
                        row.labeling,
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        is_global_winner,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )
        else:
            print(
                "\t".join(
                    [
                        row.metabolite + M_suffixes[M_offset],
                        row.formulas + M_suffixes[M_offset],
                        row.inchikey + M_suffixes[M_offset],
                        "%d" % row.fdr,
                        row.labeling,
                        row.ion_type,
                        "%.1f" % row.rt_start,
                        "%.1f" % row.rt_stop,
                        "%.1f" % (row.mz_tol * 1000000.0),
                        "%.1f" % row.rt_tol,
                        "%.4f" % abs(M_mzs[M_offset]),
                        "%.4f" % abs(local_winner_mz),
                        "%.1f" % ppm,
                        local_winner_sample,
                        is_global_winner,
                        "%.2f" % local_winner_rt,
                        "%.2f" % rt_min,
                        "%.2f" % rt_max,
                        "%.2f" % rt_range,
                        "%d" % detections,
                    ]
                    + finalized
                ),
                file=output,
            )
    output.flush()


for row in rows:
    print(f"Processing {row.metabolite}: ", file=sys.stderr, end="", flush=True)
    if ("rt" not in vars(row)) and (not (("rt_start" in vars(row)) and ("rt_stop" in vars(row)))):
        print(f"--> skipping row with insufficient information!!!", file=sys.stderr, flush=True)
        continue
    before = time.time()
    if not (("rt_start" in vars(row)) and ("rt_stop" in vars(row))):
        row["rt_start"] = max(row["rt"] - (RT_WINDOW / 2), 0.0)
        row["rt_stop"] = row["rt"] + (RT_WINDOW / 2)
    if "mz_tol" not in vars(row):
        row["mz_tol"] = MZ_TOLERANCE
    row.mz_tol = row.mz_tol / 1000000.0
    if "rt_tol" not in vars(row):
        row["rt_tol"] = RT_TOLERANCE
    is_molecular_formula = True
    c_counter = 0
    n_counter = 0
    o_counter = 0
    try:
        row.mz = float(row.formula)
        is_molecular_formula = False
        if row.mz > 0.0:
            row.polarity = "+"
        else:
            row.polarity = "-"
    except ValueError:
        # (the_mz, c_counter, n_counter, o_counter) = parse_formula(row.formula)
        parsed_formula = formula.parse(row.formula, actions=formula_actions.Actions())
        the_mz = 0.0
        c_counter = 0
        n_counter = 0
        o_counter = 0
        for term in parsed_formula:
            the_mz += term["mass"]
            if term["atom"] == "C":
                c_counter += term["count"]
            if term["atom"] == "O":
                o_counter += term["count"]
            if term["atom"] == "N":
                n_counter += term["count"]
        parsed = nist_ion_descriptions.parse(row.ion_type, actions=actions.Actions())
        row.polarity = row.ion_type[-1]
        charge = 0
        z_string = parsed["z"]
        if z_string == "+":
            charge = 1
        elif z_string == "-":
            charge = -1
        else:
            if row.polarity == "+":
                charge = int(z_string[:(-1)])
            else:
                charge = (-1) * int(z_string[:(-1)])
        row.mz = (
            (the_mz * parsed["molecular_ion_count"])
            + parsed["delta"]
            - (charge * electron)
        ) / charge
    if "labeling" in unmatched:
        exics_list = [
            exics(row.polarity, row.mz, row.rt_start, row.rt_stop, row.mz_tol)
        ]
        process_entry((row, [""], [row.mz], exics_list))
    else:
        (M_suffixes, M_values) = parse_labeling_requirements(
            row.metabolite,
            row.labeling,
            row.mz,
            is_molecular_formula,
            c_counter,
            n_counter,
            o_counter,
        )
        if not M_suffixes:
            exics_list = [
                exics(row.polarity, row.mz, row.rt_start, row.rt_stop, row.mz_tol)
            ]
            process_entry((row, [""], [row.mz], exics_list))
        else:
            exics_list = list(
                map(
                    lambda m_value: exics(
                        row.polarity, m_value, row.rt_start, row.rt_stop, row.mz_tol
                    ),
                    M_values,
                )
            )
            process_entry((row, M_suffixes, M_values, exics_list))
    after = time.time()
    print(
        f"{len(exics_list)} GICs in: {after - before :.2f} seconds ({(after - before) / len(exics_list) :.2f} / GIC)...",
        file=sys.stderr, flush=True
    )

output.close()

stop_time = time.time()

print(
    f"Skeleton processed {len(rows)} entries in {stop_time - start_time :.2f} seconds.",
    file=sys.stderr, flush=True
)
