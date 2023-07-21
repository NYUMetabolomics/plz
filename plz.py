import sys
import os
import bisect
import time
import sqlite3
import collections  # for now, plz commands are not allowed to 'from collections import something_specific'
import re
import subprocess
import operator

# import dash
# import dash_cytoscape
# import webbrowser
# import threading
# import flask

from parsers import actions
from parsers import formula_actions
from parsers import nist_ion_descriptions
from parsers import formula

import columns

#  The imports above are required, at the very least, by the command scripts (and must be passed to namespace),
#  whereas the ones below are necessitated only by plz itself...

from prompt_toolkit import print_formatted_text, HTML
from prompt_toolkit.completion import Completer, Completion
from prompt_toolkit import PromptSession
if sys.platform != "linux":
    from gooey import Gooey, GooeyParser
import hashlib


def floatable(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def intable(value):
    try:
        int(value)
        return True
    except ValueError:
        return False


if "--" in sys.argv:  # No clue why gooey does this...
    sys.argv.remove("--")

launch = True
if "--ignore-gooey" in sys.argv:
    launch = False
    sys.argv.remove("--ignore-gooey")
if len(sys.argv) >= 2:
    launch = False

if sys.platform == "linux" and len(sys.argv) == 1:
    print("To simply server setup, plz does not currently provide a GUI on linux...", file=sys.stderr)
    sys.exit(-1)

frozen = "not"
if getattr(sys, "frozen", False):
    # we are running in a bundle
    frozen = "ever so"
    bundle_dir = sys._MEIPASS
else:
    # we are running in a normal Python environment
    bundle_dir = os.path.dirname(os.path.abspath(__file__))

command_list = []
if len(sys.argv) == 2:
    the_arg = sys.argv[1]
    # for the time being, .plz files are not being developed...
    # if the_arg.endswith(".plz"):
    #     f = open(the_arg)
    #     for l in f:
    #         words = l.strip().split()
    #         command_list.append(words)
    #     f.close()
    # else:
    #     os.chdir(the_arg)
    os.chdir(the_arg)
if len(sys.argv) > 2:
    words = []
    target_dir = None
    for val in sys.argv[1:]:
        if "\\" in val:
            val = val.replace("\\", "/")
        if "/" in val:
            components = val.split("/")
            path_part = "/".join(components[:-1])
            val = components[-1]
            if not target_dir:
                target_dir = path_part
        words.append(val)
    # The first file argument dictates the script's working directory!
    if target_dir:
        os.chdir(target_dir)
    else:
        print(
            "current plz command does not specify an implicit working directory!!!",
            file=sys.stderr,
        )
    command_list.append(words)

magic_words = ["exit", "ls"]
templates = [
    template
    for template in os.listdir(bundle_dir + "/" + "commands")
    if template.endswith(".py")
]
commands = {}
for template in templates:
    code = open(bundle_dir + "/" + "commands/" + template)
    body = ""
    args = []
    imports = []
    for line in code:
        if line.startswith("#"):
            vals = line.strip().split(":")
            if vals[0] == "# output":
                continue
            the_val = vals[1].strip()
            sub_vals = the_val.split()
            #
            # TODO: factor out explicit command types.
            # The plz command template language has explicit types that cannot have defaults
            # and implicit ones that are implied by their default values. The code needs to
            # be restructured such that there is one source of truth for valid explicit types.
            # In the long term this will require a new grammar...
            #
            if floatable(the_val) or intable(the_val):
                sub_vals = []
            if the_val.startswith("."):  # used to be: the_val in ["sqlite3", "features", "quantified"]:
                sub_vals = []
            if len(sub_vals):
                the_val = sub_vals[0]
            args.append((vals[0][2:], the_val, sub_vals))
            continue
        if line.startswith("import "):
            imports.append(line.strip().split()[1])
            continue
        break
    for line in code:
        body += line
    commands[template[:(-3)]] = (args, body, imports)
for mw in magic_words:
    commands[mw] = ([], "", [])


if sys.platform != "linux":
    @Gooey(
        program_name="PLZ",
        required_cols=1,
        menu=[
            {
                "name": "File",
                "items": [
                    {
                        "type": "AboutDialog",
                        "menuTitle": "About",
                        "name": "plz",
                        "description": "Python-Based Environment for Reproducible Research Pipelines",
                        "version": "0.4.0",
                        "copyright": "Copyright © 2021 NYU Langone Health",
                        "website": "https://metabolomics.org/plz",
                        "developer": "Manor Askenazi",
                        "license": "TBD",
                    },
                    {
                        "type": "Link",
                        "menuTitle": "Visit Our Lab!",
                        "url": "https://med.nyu.edu/research/scientific-cores-shared-resources/metabolomics-laboratory/leadership",
                    },
                ],
            },
            {
                "name": "Help",
                "items": [
                    {
                        "type": "Link",
                        "menuTitle": "Documentation",
                        "url": "https://github.com/NYUMetabolomics/metorg/wiki/PLZ-Help-Page",
                    }
                ],
            },
        ],
        image_dir=bundle_dir + "/images",
    )
    def main():
        parser = GooeyParser(description="I Can Haz PipeLineZ?")
        subs = parser.add_subparsers(help="commands", dest="command")
        arg_counter = 0
        for command in commands.keys():
            if command in magic_words:
                continue
            (args, _, _) = commands[command]
            a_parser = subs.add_parser(command)
            # search_group = a_parser.add_argument_group(command)
            for (arg, argtype, constraints) in args:
                if intable(
                    argtype
                ):  # Note that the order matters (e.g. 2 is technically floatable)
                    a_parser.add_argument(arg, type=int, default=int(argtype))
                elif floatable(argtype):
                    a_parser.add_argument(arg, type=float, default=float(argtype))
                elif argtype.startswith("."):
                    a_parser.add_argument(arg, type=str, widget="FileChooser", gooey_options={
                        'wildcard': f"{argtype[1:]} files|*{argtype}"
                    })
                # elif argtype == "sqlite3":
                #     a_parser.add_argument(arg, type=str, widget="FileChooser")
                # elif argtype == "features":
                #     a_parser.add_argument(arg, type=str, widget="FileChooser")
                # elif argtype == "quantified":
                #     a_parser.add_argument(arg, type=str, widget="FileChooser")
                elif constraints:
                    a_parser.add_argument(
                        arg,
                        widget="Dropdown",
                        choices=constraints,
                        default=constraints[0],
                        gooey_options={
                            "validator": {
                                "test": 'user_input != "Select Option"',
                                "message": "Choose a value from the list",
                            }
                        },
                    )
                else:
                    a_parser.add_argument(arg, type=str)
                arg_counter += 1

        vals = parser.parse_args()
        args = vars(vals)
        return
        words = [] * arg_counter

        for arg in args:
            val = args[arg]
            if "\\" in val:
                val = val.replace("\\", "/")
            if "/" in val:
                val = val.split("/")[-1]
            words.append(val)
        execute_command(words)


keep = None


class MyCustomCompleter(Completer):
    def get_completions(self, document, complete_event):
        global keep
        keep = document
        if not document.current_line:
            for command in commands.keys():
                if len(commands[command][0]) > 0:
                    yield Completion(
                        command,
                        start_position=0,
                        display=HTML(f"<blue>{command}</blue>"),
                    )
                else:
                    yield Completion(
                        command, start_position=0, display=HTML(f"<red>{command}</red>")
                    )
        else:
            words = document.current_line.strip().split(" ")
            if len(words) == 1 and not document.current_line.endswith(" "):
                for command in commands.keys():
                    if command.startswith(words[0]):
                        if len(commands[command][0]) > 0:
                            yield Completion(
                                command[len(words[0]) :],
                                start_position=0,
                                display=HTML(f"<blue>{command}</blue>"),
                            )
                        else:
                            yield Completion(
                                command[len(words[0]) :],
                                start_position=0,
                                display=HTML(f"<red>{command}</red>"),
                            )
            if (
                len(words) == 1
                and document.current_line.endswith(" ")
                and words[0] not in commands.keys()
            ):
                return
            if (
                len(words) >= 1
                and document.current_line.endswith(" ")
                and words[0] in commands.keys()
            ):
                (args, command_code, c_imps) = commands[words[0]]
                current_arg = len(words)
                if current_arg > len(args):
                    return
                argname = args[current_arg - 1][0]
                argtype = args[current_arg - 1][1]
                constraints = args[current_arg - 1][2]
                if intable(argtype) or floatable(argtype):
                    yield Completion(
                        argtype,
                        start_position=0,
                        display=HTML(
                            f"<purple>{argtype}</purple>:<green>{argname}</green>"
                        ),
                    )
                if argtype.startswith("."):
                    sfs = [
                        suffix_file
                        for suffix_file in os.listdir(".")
                        if suffix_file.endswith(argtype)
                    ]
                    if not sfs:
                        print(f"NO {argtype} files found!")
                    for sf in sfs:
                        yield Completion(
                            sf,
                            start_position=0,
                            display=HTML(
                                f"<purple>{sf[:-8]}</purple>:<green>{argname}</green><brown>({argtype})</brown>"
                            ),
                        )
                # if argtype == "features":
                #     ffs = [
                #         feature_file
                #         for feature_file in os.listdir(".")
                #         if feature_file.endswith(".features")
                #     ]
                #     if not ffs:
                #         print("NO FEATURE FILES!")
                #     for ff in ffs:
                #         yield Completion(
                #             ff,
                #             start_position=0,
                #             display=HTML(
                #                 f"<purple>{ff[:-8]}</purple>:<green>{argname}</green><brown>({argtype})</brown>"
                #             ),
                #         )
                if constraints:
                    for option in constraints:
                        yield Completion(
                            option,
                            start_position=0,
                            display=HTML(
                                f"<purple>{option}</purple>:<green>{argname}</green>"
                            ),
                        )


def execute_command(words, soft_exit=False):
    cmd = words[0]
    cmd_complete_code = ""

    if cmd not in commands:
        print(f"Invalid command: {cmd}", file=sys.stderr)
        if soft_exit:
            return
        else:
            sys.exit(-1)
    (cmd_args, cmd_code, cmd_imports) = commands[cmd]
    import_code = ""
    for imp in cmd_imports:
        import_code += f"import {imp}\n"

    if len(cmd_args) + 1 != len(words):
        print(f"Missing arguments for command {cmd}: {len(cmd_args)} expected, but given only {len(words) - 1}")
        if soft_exit:
            return
        else:
            sys.exit(-1)
    for i, (arg, argtype, constraints) in enumerate(cmd_args):
        if intable(argtype):
            if not intable(words[i + 1]):
                print(f'"{words[i + 1]}" is not a valid integer value for the "{arg}" argument', file=sys.stderr)
                if soft_exit:
                    return
                else:
                    sys.exit(-1)
            else:
                cmd_complete_code += f'{arg} = {words[i + 1]}\n'
        elif floatable(argtype):
            if not floatable(words[i + 1]):
                print(f'"{words[i + 1]}" is not a valid floating point value for the "{arg}" argument', file=sys.stderr)
                if soft_exit:
                    return
                else:
                    sys.exit(-1)
            else:
                cmd_complete_code += f'{arg} = {words[i + 1]}\n'
        elif constraints:
            if words[i + 1] not in constraints:
                print(f'"{words[i + 1]}" not in the valid options for "{arg}" ({" ".join(constraints)})', file=sys.stderr)
                if soft_exit:
                    return
                else:
                    sys.exit(-1)
            else:
                cmd_complete_code += f'{arg} = "{words[i + 1]}"\n'
        else:
            # Since we have no 'free' string type, this means the only option left is a file type...
            # so let us also calculate a hash for it and embed that into the command.
            cmd_complete_code += f'{arg} = "{words[i + 1]}"\n'
            cmd_complete_code += f'# md5 hash of {words[i + 1]} at time of execution was {hashlib.md5(open(words[i + 1], "rb").read()).hexdigest()}\n'
    cmd_complete_code += cmd_code
    hashed_code = import_code + cmd_complete_code
    hash_object = hashlib.md5(hashed_code.encode("utf-8"))
    hash_val = hash_object.hexdigest()
    new_cmd_file_name = words[0] + "_" + hash_val + ".py"
    if words[0] != "dagger" and os.path.exists(new_cmd_file_name + ".done"):
        print("Already executed...", file=sys.stderr)
        return
    new_cmd_file = open(new_cmd_file_name, "w")
    print(hashed_code, file=new_cmd_file)
    new_cmd_file.close()
    cmd_env = {}
    cmd_env["__file__"] = new_cmd_file_name
    if "sys" in cmd_imports:
        sys_maker = collections.namedtuple("sysobj", "stderr, exit")  # TODO: Why so restrictive?
        cmd_env["sys"] = sys_maker(sys.stderr, sys.exit)
    if "os" in cmd_imports:
        cmd_env["os"] = os
    if "time" in cmd_imports:
        cmd_env["time"] = time
    if "bisect" in cmd_imports:
        cmd_env["bisect"] = bisect
    if "sqlite3" in cmd_imports:
        cmd_env["sqlite3"] = sqlite3
    if "collections" in cmd_imports:
        cmd_env["collections"] = collections
    if "re" in cmd_imports:
        cmd_env["re"] = re
    if "subprocess" in cmd_imports:
        cmd_env["subprocess"] = subprocess
    if "operator" in cmd_imports:
        cmd_env["operator"] = operator
    if "actions" in cmd_imports:
        cmd_env["actions"] = actions
    if "formula_actions" in cmd_imports:
        cmd_env["formula_actions"] = formula_actions
    if "nist_ion_descriptions" in cmd_imports:
        cmd_env["nist_ion_descriptions"] = nist_ion_descriptions
    if "formula" in cmd_imports:
        cmd_env["formula"] = formula
    if "columns" in cmd_imports:
        cmd_env["columns"] = columns
    # if "dash" in cmd_imports:
    #     cmd_env["dash"] = dash
    # if "dash_cytoscape" in cmd_imports:
    #     cmd_env["dash_cytoscape"] = dash_cytoscape
    # if "webbrowser" in cmd_imports:
    #     cmd_env["webbrowser"] = webbrowser
    # if "threading" in cmd_imports:
    #     cmd_env["threading"] = threading
    # if "flask" in cmd_imports:
    #     cmd_env["flask"] = flask
    if "mspepsearch" in cmd_imports:
        mspepsearch_maker = collections.namedtuple("mspepsearchobj", "path platform NIST Metlin LipidBLAST decoy_NIST decoy_Metlin decoy_LipidBLAST")
        if sys.platform == "win32":
            cmd_env["mspepsearch"] = mspepsearch_maker(bundle_dir + "/mspepsearch/MSPepSearch64.exe", sys.platform, "C:\\NIST14\\MSSEARCH\\nist_msms", "C:\\NIST14\\MSSEARCH\\METLIN_EXPERIMENTAL", "C:\\NIST14\\MSSEARCH\\LipidBlast_MBX", "C:\\NIST14\\MSSEARCH\\DECOY_nist_msms", "C:\\NIST14\\MSSEARCH\\DECOY_METLIN_EXPERIMENTAL", "C:\\NIST14\\MSSEARCH\\DECOY_LipidBlast_MBX")
        else:
            cmd_env["mspepsearch"] = mspepsearch_maker("/data/mspepsearch/mspepsearch", sys.platform, "/data/nist_msms", "/data/METLIN_EXPERIMENTAL", "/data/LipidBlast_MBX", "/data/DECOY_nist_msms", "/data/DECOY_METLIN_EXPERIMENTAL", "/data/DECOY_LipidBlast_MBX")
    if "lib2nist" in cmd_imports:
        lib2nist_maker = collections.namedtuple("lib2nistobj", "path platform")
        if sys.platform == "win32":
            cmd_env["lib2nist"] = lib2nist_maker(bundle_dir + "/lib2nist/lib2nist64.exe", sys.platform)
        else:
            cmd_env["lib2nist"] = lib2nist_maker(None, sys.platform)

    exec(cmd_complete_code, cmd_env)
    with open(new_cmd_file_name + ".done", "w") as finito:
        print("Finished!", file=finito)


if __name__ == "__main__":
    if launch:
        main()
    else:
        if command_list:
            for words in command_list:
                execute_command(words)
        else:
            session = PromptSession()
            while True:
                text = session.prompt(
                    HTML("<green>plz></green> "), completer=MyCustomCompleter()
                )
                words = text.split()
                cmd = words[0]
                if cmd == "exit":
                    print_formatted_text(
                        HTML("<yellow><u>You're Welcome!</u></yellow>")
                    )
                    sys.exit()
                if cmd == "ls":
                    for fname in os.listdir("."):
                        print(fname)
                    continue
                if cmd not in commands.keys():
                    print("Invalid Command!!!")
                    continue
                execute_command(words, True)
