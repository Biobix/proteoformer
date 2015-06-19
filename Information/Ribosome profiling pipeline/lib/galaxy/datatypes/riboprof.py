"""
    Ribo profiling tool datatypes
    Biobix - Ghent University
    for Ribo_prof_pipeline
    """

from galaxy import eggs

import pkg_resources
pkg_resources.require( "bx-python" )

import logging, os, sys, time, sets, tempfile, shutil
import data
from galaxy import util
from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes.binary import Binary

log = logging.getLogger(__name__)

## Ribo prof pipeline Classes

class TisOverview(Tabular):
    file_ext = 'tis'
    MetadataElement( name="args", default=[], desc="Tis Args", readonly=True, visible=True, no_value=[] )
    def __init__(self, **kwd):
        Tabular.__init__( self, **kwd )
        self.column_names = ['ID','local_max','min_count_aTIS','R_aTis','min_count_5UTR','R_5UTR','min_count_CDS','R_CDS','min_count_3UTR','R_3UTR','min_count_no_trans','R_no_trans']
        self.columns = 12
    
    def set_meta( self, dataset, overwrite = True, skip = None, max_data_lines = None, **kwd ):
        Tabular.set_meta(self, dataset, overwrite, skip, max_data_lines)
        tis_args = set()
        try:
            fh = open( dataset.file_name )
            
            for line in fh:
                fields = line.strip().split('\t')
                try:
                    tis_args.add(fields[0])
                except IndexError:
                    pass
            dataset.metadata.args = []
            dataset.metadata.args += tis_args
        
        
        finally:
            fh.close()

class SqliteDb( Binary ):
    file_ext = 'sqlitedb'
    is_binary = True
#allow_datatype_change = False



