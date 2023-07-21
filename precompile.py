import os
import sys
import shutil


def clean_dir_path(location):
    return os.path.dirname(os.path.realpath(location)).replace("\\", "/")


base_path = clean_dir_path(__file__)
conda_path = clean_dir_path(sys.executable)

print("plz path:", base_path)
print("conda path", conda_path)
shutil.copyfile("images/favicon.ico", conda_path + "/Lib/site-packages/dash/favicon.ico")
with open("plz.tspec") as f:
    with open("plz.spec", 'w') as o:
        o.write(f.read().replace("TSPECBASE", base_path).replace("TSPECCONDA", conda_path))
