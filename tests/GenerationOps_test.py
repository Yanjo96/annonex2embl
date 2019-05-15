#!/usr/bin/env python
'''
Unit Tests for the classes of the module `GenerationOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import Bio # Do not remove; important for assertIsInstance
import unittest

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

import GenerationOps as GnOps

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio import SeqFeature

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2017.02.01.1100'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()

###########
# CLASSES #
###########

class GenerateFeatLocTestCases(unittest.TestCase):
    ''' Tests for class `GenerateFeatLoc` '''

    def test_GenerateFeatLoc__make_location__1(self):
        ''' Test to evaluate function `make_location` of class `GenerateFeatLoc`.
            This test evaluates the case of a continuous range. '''
        charset_range = range(1,8)
        out = GnOps.GenerateFeatLoc().make_location(charset_range)
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation)

    def test_GenerateFeatLoc__make_location__2(self):
        ''' Test to evaluate function `make_location` of class `GenerateFeatLoc`.
            This test evaluates the case of a discontinuous range, resulting in
            a compound location. '''
        charset_range = [1,2,3,7,8]
        out = GnOps.GenerateFeatLoc().make_location(charset_range)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation)

    def test_GenerateFeatLoc__make_location__3(self):
        ''' Test to evaluate function `make_location` of class `GenerateFeatLoc`.
            This test evaluates the case of a discontinuous range that is
            separated only by a nucleotide, resulting in a compound location.
            This test evaluates if the function correctly generates a compound
            feature location if only a single base separates them. '''
        charset_range = [1,2,3,5,6]
        out = GnOps.GenerateFeatLoc().make_location(charset_range)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation) # CompoundLocation
        self.assertEqual(len(out.parts), 2)

    def test_GenerateFeatLoc__make_start_fuzzy__1(self):
        ''' Test to evaluate function `make_start_fuzzy` of class `GenerateFeatLoc`.
            This test evaluates the case where FeatureLocations are made fuzzy. '''
        from Bio import SeqFeature
        start_pos = SeqFeature.ExactPosition(5)
        end_pos = SeqFeature.ExactPosition(9)
        location_object = SeqFeature.FeatureLocation(start_pos, end_pos)
        out = GnOps.GenerateFeatLoc().make_start_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation) # FeatureLocation
        self.assertIsInstance(out.start, Bio.SeqFeature.BeforePosition) # Fuzzy Start

    def test_GenerateFeatLoc__make_start_fuzzy__2(self):
        ''' Test to evaluate function `make_start_fuzzy` of class `GenerateFeatLoc`.
            This test evaluates if start FeatureLocations are made fuzzy with a
            discontinuous location. '''
        from Bio import SeqFeature
        charset_range = [1,2,3,7,8]
        location_object = GnOps.GenerateFeatLoc().make_location(charset_range)
        out = GnOps.GenerateFeatLoc().make_start_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation) # CompoundLocation
        self.assertIsInstance(out.parts[0].start, Bio.SeqFeature.BeforePosition) # Fuzzy Start

    def test_GenerateFeatLoc__make_start_fuzzy__3(self):
        ''' Test to evaluate function `make_start_fuzzy` of class `GenerateFeatLoc`.
            This test evaluates if end FeatureLocations are made fuzzy.
            See AfterPosition. '''
        from Bio import SeqFeature
        start_pos = SeqFeature.ExactPosition(5)
        end_pos = SeqFeature.ExactPosition(9)
        location_object = SeqFeature.FeatureLocation(start_pos, end_pos)
        out = GnOps.GenerateFeatLoc().make_end_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation) # FeatureLocation
        self.assertIsInstance(out.end, Bio.SeqFeature.AfterPosition) # Fuzzy End

    def test_GenerateFeatLoc__make_start_fuzzy__4(self):
        ''' Test to evaluate function `make_start_fuzzy` of class `GenerateFeatLoc`.
            This test evaluates if end FeatureLocations are made fuzzy with a
            discontinuous location.
            See BeforePosition. '''
        from Bio import SeqFeature
        charset_range = [1,2,3,7,8]
        location_object = GnOps.GenerateFeatLoc().make_location(charset_range)
        out = GnOps.GenerateFeatLoc().make_start_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation) # CompoundLocation
        self.assertIsInstance(out.parts[0].start, Bio.SeqFeature.BeforePosition) # Fuzzy Start


class GenerateSeqFeatureTestCases(unittest.TestCase):
    ''' Tests for class `GenerateSeqFeature` '''

    def test_GenerateSeqFeature__source_feat__1(self):
        ''' Test to evaluate function `source_feat` of class `GenerateSeqFeature`.
            This test evaluates the correct generation of the SeqFeature
            `source`. A translation table is to be included. '''
        full_len = 509
        quals = {'isolate': 'taxon_B', 'country': 'Ecuador'}
        charset_names = ['foo_gene', 'foo_CDS']
        out = GnOps.GenerateSeqFeature().source_feat(full_len,
            quals, charset_names)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)

    def test_GenerateSeqFeature__source_feat__2(self):
        ''' Test to evaluate function `source_feat` of class `GenerateSeqFeature`.
            This test evaluates the correct generation of the SeqFeature
            `source`. A translation table is NOT to be included. '''
        full_len = 509
        quals = {'isolate': 'taxon_B', 'country': 'Ecuador'}
        charset_names = ['foo_intron', 'foo_IGS']
        out = GnOps.GenerateSeqFeature().source_feat(full_len,
            quals, charset_names)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)


    def test_GenerateSeqFeature__regular_feat__1(self):
        ''' Test to evaluate function `regular_feat` of class `GenerateSeqFeature`.
            This test evaluates the correct generation of a regular, non-coding
            SeqFeature. '''
        feature_name = 'psbI'
        feature_type = 'intron'
        feature_orient = 'forw'
        transl_table = 11
        feature_seq = 'ACGTACGTACGTACGT'
        charset_range = [2,3,4,5]
        feature_loc = GnOps.GenerateFeatLoc().make_location(charset_range)
        out = GnOps.GenerateSeqFeature().regular_feat(feature_name,
            feature_type, feature_orient, feature_loc, transl_table, feature_seq)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
