# unrefined: .features
# filter_by: score mcrl_score dot revdot prob percentile
# equivalence_by: inchikey inchik name mzrt 
# FDR: True False
# MAX_FDR: 100
# IF_NO_FDR_WHAT_MIN_SCORE: 900
# RT_TOL: 2.0
# PPM: 20.0
# MatchPolarity: True False
# Keep_Decoy_Hits: False True
# output: .features
import sys
import time
import operator
import actions
import nist_ion_descriptions
import columns

start_time = time.time()

__version__ = "0.2"

HeaderName = {"dot": "Dot", "revdot": "RevDot", "mcrl_score": "MCRL_Score", "score": "Score", "prob": "Prob", "percentile": "Percentile"}

(rows, unmatched, unexpected) = columns.loader(
    unrefined,
    {
        "Metabolite": {
            "field": "metabolite",
            "constructor": str,
            "required": True,
            "description": "Metabolite name.",
        },
        "InChIKey": {
            "field": "inchikey",
            "constructor": str,
            "required": False,
            "description": "InChIKey.",
        },
        "Ion Type": {
            "field": "ion_type",
            "constructor": str,
            "required": True,
            "description": "Description of ion type using the grammar described in parsers/nist_ion_descriptions.py.",
        },
        "MCRL_Score": {
            "field": "mcrl_score",
            "constructor": int,
            "required": False,
            "description": "MCRL-defined score associated with identification.",
        },
        "Score": {
            "field": "score",
            "constructor": int,
            "required": False,
            "description": "Score associated with identification.",
        },
        "Dot": {
            "field": "dot",
            "constructor": int,
            "required": False,
            "description": "Dot product based score associated with identification.",
        },
        "RevDot": {
            "field": "revdot",
            "constructor": int,
            "required": False,
            "description": "Reverse dot product based score associated with identification.",
        },
        "Prob": {
            "field": "prob",
            "constructor": float,
            "required": False,
            "description": "Probability associated with identification.",
        },
        "Source": {
            "field": "source",
            "constructor": str,
            "required": False,
            "description": "Scan identifier for the scan that led to the identification.",
        },
        "Formula": {
            "field": "formula",
            "constructor": str,
            "required": True,
            "description": "Either the Chemical Formula or the m/z of the feature identified."
        },
        "Percentile": {
            "field": "percentile",
            "constructor": float,
            "required": False,
            "description": "Highest intensity percentile of feature observation across the study/batch (higher is more intense)."

        }, 
        "RT (min)": {
            "field": "rt",
            "constructor": float,
            "required": True,
            "description": "RT (in minutes).",
        },
    },
)


PPM = PPM / 1000000.0

out = open(__file__[:-3] + ".features", "w")

if "inchikey" in unmatched:
    out_headers = "\t".join(["Metabolite", "Ion Type", "Formula", "RT (min)", HeaderName[filter_by]])
else:
    out_headers = "\t".join(["Metabolite", "InChIKey", "Ion Type", "Formula", "RT (min)", HeaderName[filter_by]])

rows.sort(key=operator.attrgetter(filter_by), reverse=True)

if equivalence_by == "mzrt":
    for row in rows:
        row.formula = float(row.formula)

ids = []
for row in rows:
    if MatchPolarity == "True":
        if row.source.split(" ")[0][-1] != row.ion_type[-1]:
            continue

    if row[filter_by] < IF_NO_FDR_WHAT_MIN_SCORE and FDR == "False":
        continue

    if row.ion_type:
        try:
            parsed = nist_ion_descriptions.parse(row.ion_type, actions=actions.Actions())
        except nist_ion_descriptions.ParseError:
            continue

    found_equivalent = False
    for existing in ids:
        if equivalence_by == "inchikey":  # existing entries without an inchikey are considered distinct...
            if existing.inchikey:
                if (row.inchikey == existing.inchikey) and (row.ion_type == existing.ion_type):
                    if abs(row.rt - existing.rt) <= RT_TOL:
                        found_equivalent = True
                        break
            else:  # matching can still occur by name when no inchikey is present in both entries...
                if (not existing.inchikey) and (row.metabolite == existing.metabolite) and (row.ion_type == existing.ion_type):
                    if abs(row.rt - existing.rt) <= RT_TOL:
                        found_equivalent = True
                        break                    
        if equivalence_by == "inchik":  # existing entries without an inchikey are considered distinct...
            if existing.inchikey:
                if (row.inchikey[:14] == existing.inchikey[:14]) and (row.ion_type == existing.ion_type):
                    if abs(row.rt - existing.rt) <= RT_TOL:
                        found_equivalent = True
                        break
            else:  # matching can still occur by name when no inchikey is present in both entries...
                if (not existing.inchikey) and (row.metabolite == existing.metabolite) and (row.ion_type == existing.ion_type):
                    if abs(row.rt - existing.rt) <= RT_TOL:
                        found_equivalent = True
                        break                    
        if equivalence_by == "name":
            if (row.metabolite == existing.metabolite) and (row.ion_type == existing.ion_type):
                if abs(row.rt - existing.rt) <= RT_TOL:
                    found_equivalent = True
                    break
        if equivalence_by == "mzrt":
            if abs((existing.formula - row.formula)/existing.formula) <= PPM:
                if abs(row.rt - existing.rt) <= RT_TOL:
                    found_equivalent = True
                    break
    if not found_equivalent:
        ids.append(row)


if equivalence_by == "mzrt":
    for row in rows:
        row.formula = f"{row.formula:0.4f}"

extract_score = operator.attrgetter(filter_by)

if FDR == "False":
    print(out_headers, file=out)
    reported_idrts = 0 
    for row in ids:
        reported_idrts += 1
        if "inchikey" in unmatched:
            print(f"{row.metabolite}\t{row.ion_type}\t{row.formula}\t{row.rt:0.2f}\t{extract_score(row)}", file=out)
        else:
            print(f"{row.metabolite}\t{row.inchikey}\t{row.ion_type}\t{row.formula}\t{row.rt:0.2f}\t{extract_score(row)}", file=out)
    out.close()
else:
    detections = []
    for row in ids:
        detections.append([row, 0])  # keep running decoy counter in situ...
    decounter = 0
    for i in range(len(detections)):
        if detections[i][0].metabolite.startswith("Decoy_"):
            decounter += 1
        detections[i][1] = decounter
    print(out_headers + "\tFDR", file=out)
    detections[-1][1] = detections[-1][1] / max(detections[-1][1], (len(detections) - detections[-1][1]))
    for i in range(len(detections)-2, 0, -1):
        detections[i][1] = min(detections[i+1][1], detections[i][1] / max(detections[i][1], (i+1) - detections[i][1]))
    reported_idrts = 0
    for entry in detections:
        if Keep_Decoy_Hits == "False" and entry[0].metabolite.startswith("Decoy_"):
            continue
        if entry[1] <= MAX_FDR / 100:  # MAX_FDR is in % now...
            reported_idrts += 1
            row = entry[0]
            if "inchikey" in unmatched:
                print(f"{row.metabolite}\t{row.ion_type}\t{row.formula}\t{row.rt:0.2f}\t{extract_score(row)}\t{100*entry[1]:.0f}", file=out)
            else:
                print(f"{row.metabolite}\t{row.inchikey}\t{row.ion_type}\t{row.formula}\t{row.rt:0.2f}\t{extract_score(row)}\t{100*entry[1]:.0f}", file=out)
stop_time = time.time()

print(f"\nreduced {len(rows)} rows to {reported_idrts} ID+RT combinations in {stop_time - start_time :.2f} seconds.", file=sys.stderr, flush=True)
