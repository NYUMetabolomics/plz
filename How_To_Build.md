# Build Environment for `plz`

Get Miniconda3 from here:

https://docs.conda.io/en/latest/miniconda.html

Install it... Note that this conda installation should be dedicated strictly for working on plz. This is because the `precompile.py` script will modify the installation. In a sense, we are doing an overkill variant of virtualenv by literally putting together a Python installation just for the purposes of building an executable version of `plz`. One way to remember that this conda is dedicated to `plz` is to name it: plzconda (during installation). It is a good idea not to add this miniconda to the path environment and not to register it as the default Python 3.8 (when running the installer).

One the installer has run, find the entry for it in the start menu (Anaconda Prompt) and from within that command shell install the following:

- pip install prompt_toolkit
- pip install gooey
- pip install pyinstaller

Remaining in that command shell use it to go to the plz directory in the metorg distribution, and run the following command:

`python precompile.py`

Finally, generate the executable using this command (the .spec file will have been generated from the preexisting .tspec file):

`pyinstaller --windowed --onefile plz.spec`

The executable will now be found in the `dist` folder.

Remember that, in order to run Dewey, the software expects the following libraries to be found exactly in these locations:

- C:\NIST14\MSSEARCH\nist_msms
- C:\NIST14\MSSEARCH\METLIN_EXPERIMENTAL
- C:\NIST14\MSSEARCH\LipidBlast_MBX
- C:\NIST14\MSSEARCH\decoy_msms_nist_msms
- C:\NIST14\MSSEARCH\decoy_msms_METLIN_EXPERIMENTAL
- C:\NIST14\MSSEARCH\decoy_msms_LipidBlast_MBX
