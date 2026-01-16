'''
REMO v1 regulatory element modules for GRCh38.

A data package providing module coordinates, metadata,
and cell ontology annotations for the human GRCh38 genome.
'''

from typing import TYPE_CHECKING, Optional

from remopy._data import modules, metadata, terms, cl_ids, tissues

if TYPE_CHECKING:
    from anndata import AnnData
    import polars as pl

__version__ = '1.0.0'
__all__ = ['modules', 'metadata', 'terms', 'cl_ids', 'tissues', 'quantify']


def quantify(
    fragments: str,
    min_fragments: int = 1000,
    modules: Optional['pl.DataFrame'] = None,
) -> 'AnnData':
    '''
    Quantify fragments into modules.
    
    See remopy._quantify.quantify for full documentation.
    Requires: pip install remopy[quantify]
    '''
    from remopy._quantify import quantify as _quantify
    return _quantify(fragments, min_fragments, modules)
