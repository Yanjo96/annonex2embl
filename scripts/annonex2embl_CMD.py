#!/usr/bin/env python2.7
'''
annonex2embl wrapper
'''

#####################
# IMPORT OPERATIONS #
#####################

import sys
import os

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

# IMPORTANT: TFL must be after "sys.path.append"
import Annonex2emblMain as AN2EMBLMain

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2019 Michael Gruenstaeudl'
__info__ = 'annonex2embl'
__version__ = '2019.05.15.1500'

#############
# DEBUGGING #
#############

import pdb
# pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

########
# TODO #
########

''' Include selection on topol of submission (linear [default] or circular) '''

############
# ARGPARSE #
############
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))

    ### REQUIRED ###
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

    parser.add_argument('-d',
                        '--descript',
                        help='text string characterizing the DNA alignment; Example: "chloroplast trnR-atpA intergenic spacer"',
                        default='[PLACEHOLDER]',
                        required=True)

    parser.add_argument('-e',
                        '--email',
                        help='Your email address; Example: "my.username@gmail.com"',
                        default='my.username@gmail.com',
                        required=True)

    parser.add_argument('-a',
                        '--authors',
                        help='Author names; Example: "Gruenstaeudl M.; LastName I."',
                        default='Gruenstaeudl M.; LastName I.',
                        required=True)

    parser.add_argument('-o',
                        '--outfile',
                        help='absolute path to outfile; outfile in EMBL format; Example: /path_to_output/test.embl',
                        default='/home/username/Desktop/test.embl',
                        required=True)

    ### OPTIONAL ###
    parser.add_argument('-ms',
                        '--manifeststudy',
                        help='Name of the study which appears in the manifest file',
                        default='',
                        required=False)

    parser.add_argument('-mn',
                        '--manifestname',
                        help='Name which appears in the manifest file',
                        required=False)

    parser.add_argument('-md',
                        '--manifestdescription',
                        help='Description for the manifest file',
                        default='',
                        required=False)

    parser.add_argument('--taxcheck',
                        help='A logical; Shall taxon names be checked against NCBI Taxonomy?',
                        default='False',
                        required=False)

    parser.add_argument('--linemask',
                        help='A logical; Shall the ID and the AC lines be masked for EntryUpload submissions?',
                        default='False',
                        required=False)

    parser.add_argument('--topol',
                        help='`circular` or `linear`.',
                        default='linear',
                        required=False)

    parser.add_argument('--taxdiv',
                        help='Any of the three letter codes specified in section 3.2 of the EMBL user manual.',
                        default='PLN',
                        required=False)

    parser.add_argument('--collabel',
                        #metavar='column specifying sequence names',
                        help='Name of column that specifies the sequence names.',
                        default='isolate',
                        required=False)

    parser.add_argument('--ttable',
                        #metavar='translation table',
                        help='Number of the translation table to translate coding regions with.'\
                        'For details, see: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi',
                        default='11',
                        required=False)

    parser.add_argument('--organelle',
                        #metavar='translation table',
                        help='Type of membrane-bound intracellular structure from which the sequence was obtained.'\
                        'For details, see: http://www.insdc.org/files/feature_table.html',
                        default='plastid',
                        required=False)

    parser.add_argument('--seqvers',
                        #metavar='sequence version',
                        help='An integer',
                        default='1',
                        required=False)

    parser.add_argument('--version',
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()


########
# MAIN #
########

    AN2EMBLMain.annonex2embl(   args.nexus,
                                args.csv,
                                args.descript,
                                args.email,
                                args.authors,
                                args.outfile,

                                args.manifeststudy,
                                args.manifestname,
                                args.manifestdescription,
                                args.taxcheck,
                                args.linemask,
                                args.topol,
                                args.taxdiv,
                                args.collabel,
                                args.ttable,
                                args.organelle,
                                args.seqvers )
