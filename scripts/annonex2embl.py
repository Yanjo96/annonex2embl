#!/usr/bin/env python2.7
'''
(a) Annotations
    Annotations are specified in the NEXUS file following format: 
    
    BEGIN SETS;
    CharSet [gene symbol]_[feature key] = [charset range];
    END;
    
    Only a single `gene symbol` and a single `feature key` must be specified per 
    line.

'''
# CURRENT DESIGN:
# Charset-definitions in .nex-file must indicate annotation type (e.g. 'CDS', 
# 'gene', 'rrna', 'trna') in their names.
# 
# INPUT: .nex, .csv
# REQUIREMENTS:
# (i) One of the columns in the csv-file must contain the sequence names 
#     also found in the nex-file.

#####################
# IMPORT OPERATIONS #
#####################

from Bio import SeqIO
#from Bio.Alphabet import generic_dna
#from Bio.Seq import Seq
from Bio import SeqFeature
from collections import OrderedDict

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

import MyExceptions as ME
import CheckingOps as CkOps
import DegappingOps as DgOps
import GenerationOps as GnOps
import ParsingOps as PrOps
import IOOps as IOOps


###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.05.03.1900'

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

###########
# CLASSES #
###########

#############
# FUNCTIONS #
#############

########
# TODO #
########
'''
TODO:
    (i) Include a function to check internet connectivity.
'''

########
# MAIN #
########

def annonex2embl(path_to_nex, path_to_csv, email_addr, path_to_outfile,
             out_format, seqname_col, transl_table):

# STEP 01: Open outfile
    outp_handle = open(path_to_outfile, 'a')

# STEP 02: Parse data from .nex-file
    try:
        charsets, alignm = IOOps.Inp().parse_nexus_file(path_to_nex)
    except ME.MyException as e:
        sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))

# STEP 03: Parse data from .csv-file
    try:
        qualifiers = IOOps.Inp().parse_csv_file(path_to_csv)
    except ME.MyException as e:
        sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))

# STEP 04: Do quality checks on input data
    try:
        CkOps.CheckCoord().quality_of_qualifiers(qualifiers, seqname_col)
    except ME.MyException as e:
        sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))

# STEP 05: Parse out feature key, obtain official gene name and gene product 
    charset_dict = {}
    for charset_name in charsets.keys():
        try:
            charset_sym, charset_type, charset_product = PrOps.ParseCharsetName(
                charset_name, email_addr).parse()
        except ME.MyException as e:
            sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))
        
        charset_dict[charset_name] = (charset_sym, charset_type,
            charset_product)

# STEP 06: Create a full SeqRecord for each sequence of the alignment.
    for seq_name in alignm.keys():

# i. Select current sequences and current qualifiers
        current_seq = alignm[seq_name]
        current_quals = [d for d in qualifiers\
            if d[seqname_col] == seq_name][0]

# ii. Generate the basic SeqRecord (i.e., without features or annotations)
        seq_record = GnOps.GenerateSeqRecord(current_seq,
            current_quals).base_record(seqname_col, charsets)

# iii. Degap the sequence while maintaing correct annotations, which has to 
#      occur before (!) the SeqFeature 'source' is generated.
#      Note: Charsets are identical across all sequences.
        degap_handle = DgOps.DegapButMaintainAnno(seq_record.seq, charsets)
        seq_record.seq, degapped_charsets = degap_handle.degap()

# iv. Generate SeqFeature 'source' and append to features list
        source_feature = GnOps.GenerateSeqFeature().source_feat(len(seq_record),
            current_quals, transl_table)
        seq_record.features.append(source_feature)

# STEP 07: Populate the feature keys with the charset information
#          Note: Each charset represents a dictionary that must be added in 
#          full to the list "SeqRecord.features"
        for charset_name, charset_range in degapped_charsets.items():
# i. Assign a gene product to a gene name
            charset_sym, charset_type, charset_product = charset_dict[charset_name]
# ii. Generate a regular SeqFeature and append to seq_record.features
            seq_feature = GnOps.GenerateSeqFeature().regular_feat(charset_sym,
                charset_type, charset_range, charset_product)
            seq_record.features.append(seq_feature)
# iii. Sort all seq_record.features except the first one (which 
#      constitutes the source feature) by the start position
        sorted_features = sorted(seq_record.features[1:],
            key=lambda x: x.location.start.position)
        seq_record.features = [seq_record.features[0]] + sorted_features

# STEP 08: Translate and check quality of translation
        for indx, feature in enumerate(seq_record.features):
            if feature.type == 'CDS' or feature.type == 'gene': # Check if feature coding region
                try:
                    feature = CkOps.CheckCoord().transl_and_quality_of_transl( \
                        seq_record, feature, transl_table)
                except ME.MyException as e:
                    print('%s annonex2embl WARNING: %s' % ('\n', e))
                    print(' Feature "%s" of sequence "%s" is not saved into '\
                        'output.' % (feature.id, seq_record.id))
                    seq_record.features.pop(indx)

# STEP 09: Write each completed record to file
        try:
            SeqIO.write(seq_record, outp_handle, out_format)
        except:
            sys.exit('%s annonex2embl ERROR: Problem with %s. Did not write to file.' % ('\n', seq_name))

# STEP 10: Close outfile
    outp_handle.close()


############
# ARGPARSE #
############

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    
    # Required
    parser.add_argument('-n',
                        '--nexus',
                        help='absolute path to infile; infile in NEXUS format; Example: /path_to_input/test.nex',
                        default='/home/username/Desktop/test.nex',
                        required=True)

    parser.add_argument('-c',
                        '--csv',
                        help='absolute path to infile; infile in CSV format; Example: /path_to_input/test.csv',
                        default='/home/username/Desktop/test.csv',
                        required=True)

    parser.add_argument('-e',
                        '--email',
                        help='Your email address',
                        default='my.username@gmail.com',
                        required=True)

    parser.add_argument('-o',
                        '--outfile',
                        help='absolute path to outfile; outfile in EMBL format; Example: /path_to_output/test.embl',
                        default='/home/username/Desktop/test.embl',
                        required=True)

    # Optional
    parser.add_argument('-f',
                        '--outformat',
                        help='Available arguments: embl, gb', 
                        default='embl',
                        required=False)

    parser.add_argument('-l',
                        '--label',
                        metavar='column specifying sequence names',
                        help='Name of column that specifies the sequence names.',
                        default='isolate',
                        required=False)

    parser.add_argument('-t',
                        '--table',
                        metavar='translation table',
                        help='Number of the translation table to translate coding regions with.'\
                        'For details, see: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi',
                        default='11',
                        required=False)

    parser.add_argument('--version', 
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

# Include selection on topology of submission (linear [default] or circular)

########
# MAIN #
########

    annonex2embl(args.nexus, args.csv, args.email, args.outfile, args.outformat, args.label, args.table)
