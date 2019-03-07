*annonex2embl*
==============

Converts an annotated DNA multi-sequence alignment (in [NEXUS](http://wiki.christophchamp.com/index.php?title=NEXUS_file_format) format) to an EMBL flatfile for submission to [ENA](http://www.ebi.ac.uk/ena) via the [commandline submission system](https://ena-docs.readthedocs.io/en/latest/cli_05.html) of Webin.


## INSTALLATION
```
python2 setup.py install  # Installation
python2 setup.py test     # Testing
```

## FILE PREPARATION
The annotations of a NEXUS file are specified via [SETS-block](http://hydrodictyon.eeb.uconn.edu/eebedia/index.php/Phylogenetics:_NEXUS_Format), which is located beneath a DATA-block and defines sets of characters in the DNA alignment. In such a SETS-block, every gene and every exon charset must be accompanied by one CDS charset. Other charsets can be defined unaccompanied.
#### Example of a complete SETS-BLOCK
```
BEGIN SETS;
CharSet trnK_intron = 0-928 2531-2813 2849-3152;
CharSet matK_gene = 929-2530;
CharSet matK_CDS = 929-2530;    # Accompanying the charset matK_gene
CharSet trnK_exon = 2814-2848;
CharSet trnK_CDS = 2814-2848;   # Accompanying the charset trnK_exon
CharSet psbA_gene = 3153-3200;
CharSet psbA_CDS = 3153-3200;   # Accompanying the charset psbA_gene
END;
```

## USAGE
#### Example with supplied test data
```
SCRPT=$PWD/scripts/annonex2embl_CMD.py
INPUT=examples/DNA_Alignment.nex
METAD=examples/Metadata.csv
DESCR="description of alignment"
EMAIL=your_email_here@gmail.com
AUTHR="Your_name_here"

python2 $SCRPT -n $INPUT -c $METAD -o ${INP%.nex*}.embl -d $DESCR -e $EMAIL -a $AUTHR
```

## TO DO
1. Add a function to compensate the contraction of an annotation due to the identification of an internal stop codon (see line 304 in Annonex2embl.py)

> Comment: Since "CkOps.TranslCheck().transl_and_quality_of_transl()" shortens annotations to the first internal stop codon encountered, the subsequent intron or IGS needs to be extended towards 5' to compensate. This can be a general function without a priori info passed to it. The important aspect is that only the SUBSEQUENT feature (if it is an intron or an IGS!) can be extended; all other features cannot be extended and need to produce a warning.

> Pseudocode:
```
Does a gap in the annotations exist?
  If yes, is the gap followed by an intron or an IGS?
    If yes, extend the intron or the IGS towards 5' to compensate.
    If no, print a warning.
  If no, continue without action.
```

2. Add a function to read in if a charset is forward or reverse and to adjust the info in the feature table.

3. Add a function that removes the accession number from the AC line ("AC   AC0663; SV 1; ..." --> "AC   XXX; SV 1; ...") and from the ID line ("ID   AC0663;" --> "ID   XXX;")

4. Add a function that automatically generates a [manifest file](https://ena-docs.readthedocs.io/en/latest/cli_01.html#manifest-file-types)

<!---
NOT NECCESARY AT THIS POINT
#### 0. Implement improvements of argparser (scripts/annonex2embl_CMD.py)
* Currently, the "required" and "optional" parameters are not displayed currently when calling scripts/annonex2embl_CMD.py. It incorrectly says "optional parameters" for all.
* Currently, --taxcheck requires "True" of "False" as parameters; how can I use it such that only the presence of --taxcheck indicates "True", whereas its abscence indicates "False"?
#### 0. Write GUI with similar to GUI of EMBL2checklists
--->

<!---
* NO LONGER NECESSARY: 0. Add a function that (a) reads and parses a bibtex file, extracts the citation info as well as the submitter references as from that file, and write the correctly formatted string-lines into the EMBL output file during post-processing.
--->


## DEVELOPMENT
#### Testing for development
To run the unittests outside of 'python setup.py test':
```
python -m unittest discover -s annonex2embl/tests -p "*_test.py"
```

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.
