1. put all the lipid MSPs in this directory
2. temporarily move parsers in to the directory (because I don't understand relative imports in python :-)
3. `python make_lipids.py 1> ion_types.txt 2> log.txt`
4. rename lipids.neo_msp into lipids.msp
5. get lib2nist and then: `lib2nist64.exe lipids.msp` (no other arguments should be necessary)