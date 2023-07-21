import sys
import os

from parsers import actions
from parsers import formula_actions
from parsers import nist_ion_descriptions
from parsers import formula
electron = 0.00054858

#
# Based on: https://github.com/suzaku/plain_obj
#


def make_constructor(fields):
    assignments = '\n'.join([f'    self.{f} = {f}' for f in fields])
    parameter_lists = ', '.join(fields)
    source = 'def __init__(self, %s):\n%s' % (parameter_lists, assignments)
    namespace = {}
    exec(source, namespace)
    return namespace['__init__']


def make_lt(key_field):
    source = f"def __lt__(self, other):\n    return self.{key_field} < other.{key_field}"
    namespace = {}
    exec(source, namespace)
    return namespace['__lt__']


def sortable(type_name, field_names):
    if isinstance(field_names, str):
        # names separated by whitespace and/or commas
        field_names = field_names.replace(',', ' ').split()
    return type(
        type_name,
        (object,),
        {
            '__slots__': field_names,
            '__init__': make_constructor(field_names),
            '__lt__': make_lt(field_names[0])
        }
    )


# This is a global used inside clear_entry.
# Essentially it is what classic OOP calls a Class Variable.
# It is 2 blank lines away from the definition because PEP8 :)
id_entry = 0


def clear_entry(an_entry):
    global id_entry
    id_entry += 1
    if not id_entry % 65536:
        # This allows us to use the default KEEP_IDS setting in LIB2NIST and nevertheless
        # circumvent its known bug regarding IDS at multiples of 65536.
        id_entry += 1
    an_entry.id = f"{id_entry}"
    an_entry.name = ""
    an_entry.formula = ""
    an_entry.ion_type = ""
    an_entry.mz = ""
    an_entry.source = ""
    an_entry.sid = ""


Entry = sortable("Entry", "id name formula ion_type mz source sid")


ion_types_seen = set()  # evil global... used by parse_file, accumulates implicitly.


def parse_file(f, sourcename, output):
    entry = Entry("", "", "", "", "", "", "")
    clear_entry(entry)  # Required to set entry.id=1
    entry.source = sourcename
    entrylines = []
    for line in f:
        line = line.strip()
        assert not line.startswith("PrecursorMZ:")  # We assume there was no explicit attempt to communicate the precursor m/z
        if line.startswith("ID:"):  # We want to _replace_ the ID entry...
            entry.sid = line[4:]
        else:
            entrylines.append(line)
        if line.startswith("Name:"):
            entry.name = line[6:]
        if line.startswith("Comment:"):
            fields = line[9:].split(";")
            (parent, mz) = fields[0].split()
            assert parent.startswith("Parent=")
            assert mz.startswith("Mz_exact=")
            parent = parent[7:]
            mz = mz[9:]
            assert parent == mz
            entry.mz = mz
            for field in fields:
                if not field:
                    continue
                field = field.strip()
                # print(f"{field}****")
                if field.startswith("["):
                    ion_types_seen.add(field)
                    if field == "[M-2H](2-)":
                        field = "[M-2H]2-"
                    if field == "[M-Ac-H]-":
                        field = "[M+C2H4O2-H]-"
                    entry.ion_type = field
                    continue
                likely_formula = True
                for char in field:
                    if char not in "CHNOPS0123456789":
                        # print(f"invalid char: {char}")
                        likely_formula = False
                        break
                if likely_formula:
                    entry.formula = field
                    continue
            if not entry.formula:
                print("Warning!!! None Formula Entry:", file=sys.stderr)
                print(line, file=sys.stderr)
                print("formula:", entry.formula, file=sys.stderr)
                print("ion_type:", entry.ion_type, file=sys.stderr)
                print("name:", entry.name, file=sys.stderr)
                print("mz:", entry.mz, file=sys.stderr)
            if not entry.ion_type:
                print("Warning!!! None Ion_Type Entry:", file=sys.stderr)
                print(line, file=sys.stderr)
                print("formula:", entry.formula, file=sys.stderr)
                print("ion_type:", entry.ion_type, file=sys.stderr)
                print("name:", entry.name, file=sys.stderr)
                print("mz:", entry.mz, file=sys.stderr)
                sys.exit(-1)
        if line.startswith("Num peaks:"):
            comment_line_was_seen = False
            mw_was_seen = False
            pnum = int(line[10:])
            for _ in range(pnum):
                peak = f.readline().strip()
                if peak.startswith("0.00000"):
                    print("Zero peak observed", file=sys.stderr)
                    pnum -= 1
                    continue
                entrylines.append(peak)
            print(entrylines[0], file=output)
            if entry.ion_type:
                print(f"Synon: $:03{entry.ion_type}", file=output)
            if entry.formula:
                print(f"Formula: {entry.formula}", file=output)
            else:
                print(f"Formula: {entry.ion_type[-1]}{entry.mz}", file=output)
            prev_line = entrylines[0]
            for remaining_line in entrylines[1:]:
                if remaining_line.startswith("Num peaks:"):
                    print(f"Num peaks: {pnum}", file=output)
                    prev_line = remaining_line
                    continue
                if remaining_line.startswith("MW:"):
                    prev_line = remaining_line
                    mw_was_seen = True
                    continue
                if prev_line.startswith("MW:"):
                    mw_was_seen = True
                    print(f"PrecursorMZ: {entry.mz}", file=output)
                if remaining_line.startswith("Comment:"):
                    comment_line_was_seen = True
                    print(f"ID: {entry.id}", file=output)
                    remaining_line = f"Comment: Source={entry.source}[{entry.sid}]"
                print(remaining_line, file=output)
                prev_line = remaining_line
            assert comment_line_was_seen
            assert mw_was_seen
            print("", file=output)
            f.readline()
            if entry.formula:
                parsed = nist_ion_descriptions.parse(entry.ion_type, actions=actions.Actions())
                parsed_formula = formula.parse(entry.formula, actions=formula_actions.Actions())
                the_mz = 0.0
                for term in parsed_formula:
                    the_mz += term["mass"]
                charge = 0
                z_string = parsed["z"]
                if z_string == "+":
                    charge = 1
                elif z_string == "-":
                    charge = -1
                else:
                    if z_string[-1] == "+":
                        charge = int(z_string[:(-1)])
                    else:
                        charge = (-1) * int(z_string[:(-1)])
                if charge != 0:
                    predicted_mz = ((the_mz * parsed["molecular_ion_count"]) + parsed["delta"] - (charge * electron)) / abs(charge)
                else:
                    predicted_mz = 0.0
                if abs(float(entry.mz) - predicted_mz) > 0.005:
                    print(f"new_id: {entry.id}")
                    print(f"name: {entry.name}")
                    print(f"ion_type: {entry.ion_type}")
                    print(f"formula: {entry.formula}")
                    print(f"mz: {entry.mz}")
                    print(f"source: {entry.source}")
                    print(f"source_id: {entry.sid}")
                    print(f"=====")
                    print(f"base mass: {the_mz}")
                    print(f"parsed_z: {parsed['z']}")
                    print(f"calculated charge: {charge}")
                    print(f"molecular_ion_count: {parsed['molecular_ion_count']}")
                    print(f"delta: {parsed['delta']}")
                    print(f"predicted_mz: {predicted_mz}")
                    sys.exit()
            clear_entry(entry)
            entry.source = sourcename
            entrylines = []


with open("lipids.neo_msp", 'w') as output:
    libnames = [libname for libname in os.listdir(".") if libname.endswith(".MSP")]
    for libname in libnames:
        with open(libname) as f:
            parse_file(f, libname[:-4], output)
            id_entry -= 1  # Very Dangerous -- this is just to unwind one unnecessary increment that always happens per file...
print(ion_types_seen)

# libnames = [libname for libname in os.listdir(".") if libname.endswith(".MSP") and (not libname.startswith("neo_"))]
# for libname in libnames:
#     with open("neo_" + libname, 'w') as output:
#         with open(libname) as f:
#             parse_file(f, libname[:-4], output)
