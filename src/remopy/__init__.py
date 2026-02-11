'''
REMO v1 regulatory element modules for GRCh38.

A data package providing module coordinates, metadata,
and cell ontology annotations for the human GRCh38 genome.
'''

from typing import TYPE_CHECKING

from remopy._data import modules, metadata, terms, ontology, tissues

if TYPE_CHECKING:
    from pathlib import Path
    from anndata import AnnData
    import polars as pl

__version__ = '1.0.2'
__all__ = ['modules', 'metadata', 'terms', 'ontology', 'tissues', 'quantify']


def quantify(
    fragments: 'str | Path',
    min_fragments: int = 1000,
    modules: 'pl.DataFrame | None' = None,
    verbose: bool = True,
) -> 'AnnData':
    '''
    Quantify fragments into modules.

    Counts the number of fragments overlapping each module per cell.
    Each fragment is counted once (the read support column is ignored).

    Parameters
    ----------
    fragments
        Path to fragments file (10x format: chrom, start, end, barcode, [count]).
        Can be gzipped.
    min_fragments
        Minimum total fragments per cell to include in output.
    modules
        Optional filtered modules DataFrame. If provided, only these modules
        will be used for quantification. Must have columns: chrom, start, end, REMO.
        If None, all modules are used.
    verbose
        If True, log progress messages. If False, suppress output.

    Returns
    -------
    AnnData
        AnnData object with cells as obs and modules as var.
        X contains fragment counts (sparse integer matrix).
        var contains module metadata (CREs, Bases, GC_mean, CL, etc.).

    Examples
    --------
    >>> import remopy as remo
    >>> adata = remo.quantify('fragments.tsv.gz')
    >>> adata
    AnnData object with n_obs × n_vars = 5000 × 340069

    >>> # Quantify only T cell modules
    >>> t_cell_ids = remo.terms().get('T cell', [])
    >>> t_cell_mods = remo.modules().filter(pl.col('REMO').is_in(t_cell_ids))
    >>> adata = remo.quantify('fragments.tsv.gz', modules=t_cell_mods)

    Requires
    --------
    pip install remopy[quantify]
    '''
    from remopy._quantify import quantify as _quantify
    return _quantify(fragments, min_fragments, modules, verbose)
