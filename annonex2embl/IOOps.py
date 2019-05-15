#!/usr/bin/env python
'''
Custom operations input and output processes
'''

#####################
# IMPORT OPERATIONS #
#####################

import os
import MyExceptions as ME

from csv import DictReader
from Bio.Nexus import Nexus
from StringIO import StringIO
from Bio import SeqIO

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

###########
# CLASSES #
###########

class Inp:
    ''' This class contains functions to conduct miscellaneous input
        operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def extract_fn(self, in_path):
        ''' This function splits a the path from path+filename. '''
        path, fn = os.path.split(in_path)
        return fn

    def repl_fileend(self, fn, new_end):
        ''' This function replaces the file ending (e.g. ".csv") with a
        different ending. '''
        return fn[:fn.rfind('.')] + '.' + new_end


    def parse_csv_file(self, path_to_csv):
        ''' This function parses a csv file. '''
        try:
            reader = DictReader(open(path_to_csv, 'rb'), delimiter=',',
                                quotechar='"', skipinitialspace=True)
            a_matrix = list(reader)
        except:
            raise ME.MyException('Parsing of .csv-file unsuccessful.')
        return a_matrix

    def parse_nexus_file(self, path_to_nex):
        ''' This function parses a NEXUS file. '''
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_nex)
            charsets = aln.charsets
            matrix = aln.matrix
        except Nexus.NexusError as ne:
            raise ne
        except:
            raise ME.MyException('Parsing of .nex-file unsuccessful.')
        return (charsets, matrix)


class Outp:
    ''' This class contains two functions for various output operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def write_EntryUpload(self, seq_record, outp_handle, eusubm_bool):
        ''' This function writes a seqRecord in ENA format for a submission
            via Entry Upload. Upon request (eusubm_bool), it also masks the ID and AC
            lines as requested by ENA for submissions.
        Args:
            seq_record (obj)
            outp_handle (obj)
            eusubm_bool(str)
        Returns:
            currently nothing
        Raises:
            -
        '''
        temp_handle = StringIO()
        try:
            SeqIO.write(seq_record, temp_handle, 'embl')
        except:
            raise ME.MyException('%s annonex2embl ERROR: Problem with \
            `%s`. Did not write to internal handle.' % ('\n', seq_name))
        if eusubm_bool:
            temp_handle_lines = temp_handle.getvalue().splitlines()
            if temp_handle_lines[0].split()[0] == 'ID':
                ID_line = temp_handle_lines[0]
                ID_line_parts = ID_line.split('; ')
                if len(ID_line_parts) == 7:
                    ID_line_parts = ['XXX' if ID_line_parts.index(p) in
                                     [0, 1, 3, 4, 5, 6] else p for p in ID_line_parts]
                temp_handle_lines[0] = 'ID   ' + '; '.join(ID_line_parts)
            if temp_handle_lines[2].split()[0] == 'AC':
                temp_handle_lines[2] = 'AC   XXX;'
            temp_handle_new = '\n' + '\n'.join(temp_handle_lines)
            temp_handle.truncate(0)
            temp_handle.write(temp_handle_new)
        else:
            pass

        outp_handle.write(temp_handle.getvalue())
        temp_handle.close()

        # return something?

    def create_manifest_file(self, path_to_outfile, study, name, description = ""):
        test = ''.join(path_to_outfile.split('.')[:-1]) + '.manifest'
        manifest = open(test, "w")
        manifest.write("STUDY\t" + study + "\n")
        manifest.write("NAME\t" + name + "\n")
        if(description != ""):
            manifest.write("Description\t" + description + "\n")
        manifest.close()
