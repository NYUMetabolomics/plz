# Installing `plz`

Basically, `plz.exe` is a completely self-contained executable. Just copy it onto your machine and double-click... This will _not_ install any software permanently -- it will simply run. That means you don't need to bother keeping it on your machine, just get another copy next time you need to run a pipeline on your PC!

The big exception is when running dewey, which is the library search command. For this you will need to have three libraries on your machine: NIST msms, Metlin and our very own Lipid_MBX.

Simply place the folders (found in `Official_Spectral_Libraries`) named:

- nist_msms
- METLIN_EXPERIMENTAL
- LipidBlast_MBX
- DECOY_nist_msms
- DECOY_METLIN_EXPERIMENTAL
- DECOY_LipidBlast_MBX

Into the following directory: `C:\NIST14\MSSEARCH`
