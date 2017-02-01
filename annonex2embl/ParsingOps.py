#!/usr/bin/env python
'''
Classes to parse charset names
'''

#####################
# IMPORT OPERATIONS #
#####################

import GlobalVariables as GlobVars
import MyExceptions as ME

import sys

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2017.01.31.2000'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()

###########
# CLASSES #
###########

class GetEntrezInfo:
    ''' This class contains functions to obtain gene information from gene
    symbols. '''

    def __init__(self, email_addr):
        self.email_addr = email_addr

    @staticmethod
    def _id_lookup(gene_sym, retmax=10):
        ''' An internal static function to convert a gene symbol to an Entrez ID 
            via ESearch.
        Args:
            gene_sym (str): a gene symbol; example: 'psbI'
        Returns:
            entrez_id_list (list): a list of Entrez IDs; example: ['26835430',
                            '26833718', '26833393', ...]
        Raises:
            none
        '''

#        Examples:
#            Example 1: # Default behaviour
#                >>> gene_sym = 'psbI'
#                >>> _id_lookup(gene_sym)
#                Out: ['26835430', '26833718', '26833393', ...]

        from Bio import Entrez
        if not gene_sym:
            raise ME.MyException('No gene symbol detected.')
        if '_' in gene_sym:
            raise ME.MyException('Gene symbol `%s` contains an '\
                'underscore, which is not allowed.' % (gene_sym))
        query_term = gene_sym + ' [sym]'
        try:
            esearch_records = Entrez.esearch(db='gene', term=query_term,
                retmax = retmax, retmod='xml')
        except:
            raise ME.MyException('An error occurred while retrieving '\
                'data from %s.' % ('ESearch'))
        parsed_records = Entrez.read(esearch_records)
        entrez_id_list = parsed_records['IdList']
        return entrez_id_list

    @staticmethod
    def _gene_product_lookup(entrez_id_list):
        ''' An internal static function to convert a list of Entrez IDs to a
        list of Entrez gene records via EPost and ESummary.
        Args:
            entrez_id_list (list): a list of Entrez IDs; example: ['26835430',
                                   '26833718', '26833393', ...]
        Returns:
            entrez_rec_list (list): a list of Entrez gene records
        Raises:
            none
        '''

#        Examples:
#            Example 1: # Default behaviour
#                >>> entrez_id_list = ['26835430', '26833718', '26833393']
#                >>> _record_lookup(entrez_id_list)
#                Out: ???

        from Bio import Entrez
        epost_query = Entrez.epost('gene', id=','.join(entrez_id_list))
        try:
            epost_results = Entrez.read(epost_query)
        except:
            raise ME.MyException('An error occurred while retrieving data from '\
                '%s.' % ('EPost'))
        webenv = epost_results['WebEnv']
        query_key = epost_results['QueryKey']
        try:
            esummary_records = Entrez.esummary(db='gene', webenv=webenv,
                query_key=query_key)
        except:
            raise ME.MyException('An error occurred while retrieving data from '\
                '%s.' % ('ESummary'))
        entrez_rec_list = Entrez.read(esummary_records)
        return entrez_rec_list

    @staticmethod
    def _parse_gene_products(entrez_rec_list):
        ''' An internal static function to parse out relevant information from .
        Args:
            entrez_rec_list (list): a list of Entrez gene records
        Returns:
            gene_info_list (list): a list of dictionaries
        Raises:
            none
        '''

#        Examples:
#            Example 1: # Default behaviour
#                >>> entrez_rec_list = []
#                >>> _parse_records(entrez_rec_list)

        from collections import Counter
        try:
            documentSummarySet = entrez_rec_list['DocumentSummarySet']
            docs = documentSummarySet['DocumentSummary']
        except:
            raise ME.MyException('An error occurred while parsing the '\
            'data from %s.' % ('ESummary'))
        list_gene_product = [doc['Description'] for doc in docs]
        #list_gene_symbol = [doc['NomenclatureSymbol'] for doc in docs]
        #list_gene_name = [doc['Name'] for doc in docs]

        # Avoiding that spurious first hit biases gene_product
        gene_product = Counter(list_gene_product).most_common()[0][0]
        return gene_product

    @staticmethod
    def _taxname_lookup(taxon_name, retmax=1):
        ''' An internal static function to look up a taxon name at NCBI 
            Taxonomy via ESearch.
        Args:
            taxon_name (str): a taxon name; example: 'Pyrus tamamaschjanae'
            retmax (int):     the number of maximally retained hits
        Returns:
            entrez_hitcount (int): an integer
        Raises:
            none
        '''

#        Examples:
#            Example 1: # Default behaviour
#                >>> taxon_name = 'Pyrus tamamaschjanae'
#                >>> _taxname_lookup(taxon_name)
#                Out: 0
        
        from Bio import Entrez
        if not taxon_name:
            raise ME.MyException('No taxon name detected.')
        if '_' in taxon_name:
            raise ME.MyException('Taxon name `%s` contains an underscore, '
            'which is not allowed.' % (taxon_name))
        query_term = taxon_name
        try:
            esearch_records = Entrez.esearch(db='taxonomy', term=query_term,
                retmax=retmax, retmod='xml')
        except:
            raise ME.MyException('An error occurred while retrieving data from '\
                '%s.' % ('ESearch'))
        parsed_records = Entrez.read(esearch_records)
        entrez_hitcount = parsed_records['Count']
        return entrez_hitcount


    def obtain_gene_product(self, gene_sym):
        ''' This function performs something.
        '''

#        Examples:
#            Example 1: # Default behaviour
#                >>> gene_sym = 'psbI'
#                >>> GetGeneInfo()._entrezid_lookup(gene_sym)
#                Out: ['26835430', '26833718', '26833393', ...]

        from Bio import Entrez
        Entrez.email = self.email_addr
        try:
            entrez_id_list = GetEntrezInfo._id_lookup(gene_sym)
        except ME.MyException as e:
            raise e
        try:
            entrez_rec_list = GetEntrezInfo._gene_product_lookup(entrez_id_list)
        except ME.MyException as e:
            raise e
        try:
            gene_product = GetEntrezInfo._parse_gene_products(entrez_rec_list)
        except ME.MyException as e:
            raise e
        return gene_product


    def does_taxon_exist(self, taxon_name):
        ''' This function calls _taxname_lookup and thus evaluates if a taxon exists.
        Args:
            taxon_name (str): a taxon name; example: 'Pyrus tamamaschjanae'
            retmax (int):     the number of maximally retained hits
        Returns:
            entrez_id_list (list): a list of Entrez IDs; example: ['26835430',
                            '26833718', '26833393', ...]
        Raises:
            none
        '''
        from Bio import Entrez
        Entrez.email = self.email_addr
        try:
            entrez_hitcount = GetEntrezInfo._taxname_lookup(taxon_name)
        except ME.MyException as e:
            raise e
        if entrez_hitcount == '0':
            return False
        if entrez_hitcount == '1':
            return True


class ConfirmAdjustTaxonName:
    ''' This class contains functions to confirm or adjust a sequence's
    taxon name.
    '''

    def __init__(self):
        pass
    
    def go(self, seq_record, email_addr):
        ''' This function evaluates a taxon name against NCBI taxonomy; 
            if not listed, it adjusts the taxon name and appends it
            as ecotype info.
            Args:
                seq_record (obj):   a seqRecord object
                email_addr (dict):  your email address; example: 
                                    "m.gruenstaeudl@fu-berlin.de"
            Returns:
                seq_record (obj):   a seqRecord object
            Raises:
                currently nothing
        '''
        species_name = seq_record.name
        try:
            genus_name = species_name.split()[0]
        except ME.MyException as e:
            sys.exit('%s annonex2embl ERROR: Could not locate a '\
                'whitespace between genus name and specific epithet '\
                'in taxon name of sequence `%s`.' % ('\n', seq_record.id))
        if not GetEntrezInfo(email_addr).does_taxon_exist(species_name):
            print('%s annonex2embl WARNING: Taxon name of sequence `%s` '\
                'not found in NCBI Taxonomy. Thus, adjusting taxon name.' 
                % ('\n', seq_record.id))
            if not GetEntrezInfo(email_addr).does_taxon_exist(genus_name):
                sys.exit('%s annonex2embl ERROR: Neither genus name, '\
                    'nor species name of sequence `%s` were found in '\
                    'NCBI Taxonomy.' % ('\n', seq_record.id))
            else:
                species_name_original = species_name
                species_name_new = genus_name + ' sp.'
                seq_record.name = species_name_new
                seq_record.features[0].qualifiers['organism'] = species_name_new
                seq_record.description = seq_record.description.\
                    replace(species_name_original, species_name_new)
                if 'ecotype' in seq_record.features[0].qualifiers:
                    temp = seq_record.features[0].qualifiers['ecotype']
                    seq_record.features[0].qualifiers['ecotype'] = \
                        temp + '; ' + species_name_original
                else:
                    seq_record.features[0].qualifiers['ecotype'] = \
                        species_name_original
        return seq_record


class ParseCharsetName:
    ''' This class contains functions to parse charset names. 
    Args:
        charset_name (str): a string that represents a charset name; example: 
                            "psbI_CDS"
        email_addr (dict):  your email address; example: 
                            "m.gruenstaeudl@fu-berlin.de"
    Raises:
        currently nothing
    '''

    def __init__(self, charset_name, email_addr):
        self.charset_name = charset_name
        self.email_addr = email_addr

    @staticmethod
    def _extract_charset_type(charset_name):
        ''' An internal static function to extract the charset type from a 
        string. '''
        fk_present = [fk for fk in GlobVars.nex2ena_valid_INSDC_featurekeys 
            if fk in charset_name]
        if not fk_present:
            raise ME.MyException('%s annonex2embl ERROR: No feature '\
            'key encountered in the name of charset `%s`.' % ('\n',
            charset_name))
        if len(fk_present) > 1:
            raise ME.MyException('%s annonex2embl ERROR: More than '\
            'one feature key encountered in the name of charset '\
            '`%s`.' % ('\n', charset_name))
        charset_type = fk_present[0]
        return charset_type
    
    @staticmethod
    def _extract_charset_sym(charset_name, charset_type):
        ''' An internal static function to extract the charset symbol from a 
        string. '''
        try:
            charset_sym = charset_name.strip(charset_type)
        except:
            raise ME.MyException('%s annonex2embl ERROR: No charset '\
            'symbol encountered in the name of charset `%s`.' % (
            '\n', charset_name))
        charset_sym = charset_sym.strip('_')
        charset_sym = charset_sym.rstrip('_') # Remove trailing underscores
        return charset_sym

    def parse(self):
        ''' This function parses the charset_name.
        Returns:
            tupl.   The return consists of three strings in the order 
                    "charset_sym, charset_type, charset_product"
        '''
        try:
            charset_type = ParseCharsetName._extract_charset_type(\
                self.charset_name)
        except ME.MyException as e:
            raise e
        try:
            charset_sym = ParseCharsetName._extract_charset_sym(\
                self.charset_name, charset_type)
        except ME.MyException as e:
            raise e
        entrez_handle = GetEntrezInfo(self.email_addr)
        if charset_type == 'CDS' or charset_type == 'gene':
            try:
                charset_product = entrez_handle.obtain_gene_product(\
                    charset_sym)
            except ME.MyException as e:
                raise e
        else:
            charset_product = None
        return (charset_sym, charset_type, charset_product)


#############
# FUNCTIONS #
#############

########
# MAIN #
########
