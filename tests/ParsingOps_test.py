#!/usr/bin/env python
'''
Unit Tests for the classes of the module `ParsingOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import unittest

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

import GlobalVariables as GlobVars
import MyExceptions as ME
import ParsingOps as PrOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2017.02.01.1400'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()

###########
# CLASSES #
###########


class ParseCharsetNameTestCases(unittest.TestCase):
    ''' Tests to evaluate class `ParseCharsetName` '''

    def test_ParseCharsetName__parse__1(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the basic parsing ability of the
            class. '''
        charset_name = 'psbI_CDS'
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        handle = PrOps.ParseCharsetName(charset_name, email_addr).parse()
        self.assertIsInstance(handle, tuple)
        self.assertIsInstance(handle[0], str)
        self.assertIsInstance(handle[1], str)
        self.assertTrue(handle[1] in GlobVars.nex2ena_valid_INSDC_featurekeys)

    def test_ParseCharsetName__parse__2(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name
            contains more than one charset_type. '''
        charset_name = 'psbI_rRNA_CDS' # Two feature keys present.
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()

    def test_ParseCharsetName__parse__3(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name contains
            more than one charset_sym. '''
        charset_name = 'psbI_matK_CDS' # Two gene symbols present.
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()

    def test_ParseCharsetName__parse__4(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name contains an
            unknown charset_sym. '''
        charset_name = 'xxxX_CDS'
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()

    def test_ParseCharsetName__parse__5(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name
            contains only the charset_sym, not the charset_type. '''
        charset_name = 'matK'
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()

    def test_ParseCharsetName__parse__6(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name
            contains only the charset_type, not the charset_sym. '''
        charset_name = 'CDS'
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()

    def test_ParseCharsetName__parse__7(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where the position of the
            charset_sym and the charset_type have been inadvertedly
            switched. Such a switch should not matter. '''
        charset_name = 'CDS_psbI'
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        handle = PrOps.ParseCharsetName(charset_name, email_addr).parse()
        self.assertIsInstance(handle, tuple)
        self.assertIsInstance(handle[0], str)
        self.assertIsInstance(handle[1], str)
        self.assertTrue(handle[1] in GlobVars.nex2ena_valid_INSDC_featurekeys)

class GetEntrezInfoTestCases(unittest.TestCase):
    ''' Tests to evaluate class `GetEntrezInfo` '''

    def test_GetEntrezInfo__does_taxon_exist__1(self):
        ''' This test evaluates function `does_taxon_exist` of class `GetEntrezInfo`.
            This test evaluates the case where a taxon name is used that does
            not exist on NCBI Taxonomy. '''
        #taxon_name = 'Pyrus tamamaschjanae'
        taxon_name = "qwertzuiop"
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        handle = PrOps.GetEntrezInfo(email_addr).does_taxon_exist(taxon_name)
        self.assertFalse(handle)

    def test_GetEntrezInfo__does_taxon_exist__2(self):
        ''' This test evaluates function `does_taxon_exist` of class `GetEntrezInfo`.
            This test evaluates the case where a taxon name is used that does
            exist on NCBI Taxonomy. '''
        taxon_name = 'Pyrus caucasica'
        email_addr = 'm.gruenstaeudl@fu-berlin.de'
        handle = PrOps.GetEntrezInfo(email_addr).does_taxon_exist(taxon_name)
        self.assertTrue(handle)

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
