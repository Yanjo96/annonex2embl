CHANGELOG
---------
#### Version 0.4.6
* Added function that add "/codon_start=1" in CDS feature, if start and stop position of feature is uncertain
* gene and exon features than are less than 15 nt long, the annotation would be dropped form the output
* Added a function that automatically removes all sequences that consists only of Ns (ot ?s)
* Added a function that compensate the contraction of an annotation due to the identification of an internal stop codon
* Added a function that read in if a charset is forward or reverse
* Added a function that automatically generates a [manifest file](https://ena-docs.readthedocs.io/en/latest/cli_01.html#manifest-file-types)
* Added a function that removes the accession number from the AC line and from the ID line if the user wants to
#### Version 0.4.5 (2018.05.22)
* Added function that converts missing sections of a sequence that are longer than 2 nucleotides into a "gap"-feature
#### Version 0.4.4 (2018.03.29)
* Separate the annonex2embl from the embl2checklist function.
#### Version 0.4.3 (2018.03.23)
* Improved formatting of Python code
* Checking if sequence names in NEX-file identical to sequence ids in csv-file
* Discarding charset_ranges that are empty
* Source feature 'organelle' implemented
#### Version 0.4.2 (2017.02.01)
* All qualifier values are formatted to consist of ASCII characters only.
* If a coding region is among the sequence features, the qualifier /\trans_table/ is added to the source feature.
* Implementation of customization of DE line.
#### Version 0.3 (2017.01.25)
* Rudimentary checklist output has been implemented.
* Taxon names are compared to [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) and converted to informal names if not listed there.
#### Version 0.2 (2017.01.22)
* Leading and trailing ambiguities per sequence are removed while maintaining correct annotations.
* Modifiers without content are no longer included in the output.
* Specification of taxonomic division and sequence version has been implemented.
* A submission mode that masks the ID line and the AC line has been implemented.
