'''Data loading functions.'''

from functools import cache
from importlib.resources import files
import json

import polars as pl


def _path(name: str) -> str:
    '''Get path to data file.'''
    return str(files('remopy.data').joinpath(name))


@cache
def modules() -> pl.DataFrame:
    '''
    Load REMO module genomic coordinates.
    
    Returns
    -------
    pl.DataFrame
        DataFrame with columns: chrom, start, end, REMO
    '''
    return pl.read_csv(
        _path('REMOv1_GRCh38.bed.gz'),
        separator='\t',
        has_header=False,
        new_columns=['chrom', 'start', 'end', 'REMO']
    )


@cache
def metadata() -> pl.DataFrame:
    '''
    Load REMO module metadata.
    
    Returns
    -------
    pl.DataFrame
        DataFrame with columns: REMO, CREs, Bases, Chromosome, GC_mean, CL
    '''
    return pl.read_parquet(_path('metadata.parquet'))


@cache
def terms() -> dict[str, list[str]]:
    '''
    Load Cell Ontology term to REMO module mappings.
    
    Returns
    -------
    dict
        Mapping of cell type name -> list of REMO module IDs
    '''
    with open(_path('terms.json')) as f:
        return json.load(f)


@cache
def cl_ids() -> dict[str, list[str]]:
    '''
    Load Cell Ontology ID to REMO module mappings.
    
    Returns
    -------
    dict
        Mapping of CL ID (e.g., 'CL:0000057') -> list of REMO module IDs
    '''
    with open(_path('cl_ids.json')) as f:
        return json.load(f)


@cache
def tissues() -> dict[str, list[str]]:
    '''
    Load tissue to cell type mappings.
    
    Returns
    -------
    dict
        Mapping of tissue name -> list of cell type names present in that tissue
    '''
    with open(_path('tissues.json')) as f:
        return json.load(f)
