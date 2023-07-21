import sys


class Thing:
    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)


def enum(options):
    values = options.split(" ")

    def __effectively_enums_are_strings__(element):
        if element in values:
            return str(element)
        print(f"Invalid value: {element} is not one of: {options}", sys.stderr)
        sys.exit(-1)

    return __effectively_enums_are_strings__


def loader(filename, header_description):
    for available in header_description:
        assert header_description[available]["field"] != "rest_of_row"
    f = open(filename)
    headers = f.readline().strip(" \r\n").split("\t")
    header_count = len(headers)
    to_parse = {}
    matched = {}
    matched_headers = []
    unexpected = []
    unmatched = []
    for (offset, header) in enumerate(headers):
        matched_it = False
        for available in header_description:
            if header == available:
                to_parse[offset] = header_description[available]
                matched[header_description[available]["field"]] = offset
                matched_headers.append(available)
                matched_it = True
                break
        if not matched_it:
            unexpected.append(header)
            to_parse[offset] = None
    for available in header_description:
        if available not in matched_headers:
            if header_description[available]["required"]:
                print(
                    f'Required header "{available}" was not found in the headers of {filename}!!!',
                    file=sys.stderr,
                )
                sys.exit(-1)
            else:
                unmatched.append(header_description[available]["field"])

    rows = []
    for line in f:
        new_row = Thing()
        vals = line.strip(" \r\n").split("\t")
        others = []
        for i in range(header_count):
            if to_parse[i]:
                new_row[to_parse[i]["field"]] = (to_parse[i]["constructor"])(vals[i])
            else:
                others.append(vals[i])
        if others:
            new_row["rest_of_row"] = "\t".join(others)
        rows.append(new_row)
    return (rows, unmatched, unexpected)


# Usage example:

# (rows, unmatched, unexpected) = loader(
#     input_filename,
#     {
#         "Metabolite": {
#             "field": "metabolite",
#             "constructor": str,
#             "required": True,
#             "description": "Metabolite name.",
#         },
#         "Formula": {
#             "field": "formula",
#             "constructor": str,
#             "required": True,
#             "description": "Formula or m/z value.",
#         },
#         "InChIKey": {
#             "field": "inchikey",
#             "constructor": str,
#             "required": False,
#             "description": "InChIKey.",
#         },
#         #  "Polarity (z)": {
#         #     "field": "polarity",
#         #     "constructor": enum("+ -"),
#         #     "required": True,
#         #     "description": "Polarity: (+/-)."
#         # },
#         "Labeling": {
#             "field": "labeling",
#             "constructor": str,
#             "required": False,
#             "description": "Labeling scheme assumed for the given metabolite.",
#         },
#         "Ion Type": {
#             "field": "ion_type",
#             "constructor": str,
#             "required": True,
#             "description": "Ionization type.",
#         },
#         "RT Start (min)": {
#             "field": "rt_start",
#             "constructor": float,
#             "required": True,
#             "description": "RT window start (in minutes).",
#         },
#         "RT End (min)": {
#             "field": "rt_stop",
#             "constructor": float,
#             "required": True,
#             "description": "RT window end (in minutes).",
#         },
#         "m/z Tolerance (ppm)": {
#             "field": "mz_tol",
#             "constructor": float,
#             "required": True,
#             "description": "m/z tolerance in PPM.",
#         },
#         "RT Tolerance (min)": {
#             "field": "rt_tol",
#             "constructor": float,
#             "required": True,
#             "description": "RT tolerance (in minutes).",
#         },
#     },
# )
