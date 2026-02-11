'''Fragment quantification into REMO.'''

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse
from anndata import AnnData

from remopy._data import modules as load_modules, metadata

logger = logging.getLogger(__name__)


def quantify(
    fragments: str | Path,
    min_fragments: int = 1000,
    modules: pl.DataFrame | None = None,
    verbose: bool = True,
) -> AnnData:
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
    '''
    try:
        import polars_bio as pb
    except ImportError:
        raise ImportError(
            'polars-bio is required for quantification. '
            'Install with: pip install polars-bio'
        )

    if not Path(fragments).exists():
        raise FileNotFoundError(f'Fragments file not found: {fragments}')
    if min_fragments < 1:
        raise ValueError(f'min_fragments must be >= 1, got {min_fragments}')

    if verbose:
        logger.setLevel(logging.INFO)
        if not logger.handlers:
            logger.addHandler(logging.StreamHandler())
    else:
        logger.setLevel(logging.CRITICAL)

    pb.set_option('datafusion.bio.coordinate_system_check', False)
    pb.set_option('datafusion.bio.coordinate_system_zero_based', True)

    # Load REMO modules
    mods = modules if modules is not None else load_modules()
    
    # Load fragments (only first 4 columns needed)
    frags = pl.read_csv(
        fragments,
        separator='\t',
        has_header=False,
        columns=[0, 1, 2, 3],
        new_columns=['chrom', 'start', 'end', 'barcode']
    )
    
    # Filter cells by total fragment count
    cell_counts = frags.group_by('barcode').len()
    keep = cell_counts.filter(pl.col('len') >= min_fragments).select('barcode')
    n_cells_before = cell_counts.height
    n_cells_after = keep.height
    frags = frags.join(keep, on='barcode')
    logger.info(f'Filtered to {n_cells_after}/{n_cells_before} cells with >= {min_fragments} fragments')
    
    # Find overlaps between fragments and modules
    logger.info('Finding overlaps...')
    overlaps = pb.overlap(frags, mods, suffixes=('', '_mod'), output_type='polars.DataFrame')
    del frags
    
    # Count fragments per cell per module
    logger.info('Counting fragments per module...')
    counts = overlaps.group_by(['barcode', 'REMO']).len().rename({'len': 'n'})
    del overlaps
    
    # Build sparse matrix
    logger.info('Building count matrix...')
    barcodes = counts['barcode'].unique().sort().to_list()
    remo_ids = mods['REMO'].unique().sort().to_list()
    
    bc_map = pl.DataFrame({'barcode': barcodes, 'row': range(len(barcodes))})
    mod_map = pl.DataFrame({'REMO': remo_ids, 'col': range(len(remo_ids))})
    counts = counts.join(bc_map, on='barcode').join(mod_map, on='REMO')
    
    X = sparse.csr_matrix(
        (counts['n'].to_numpy(), (counts['row'].to_numpy(), counts['col'].to_numpy())),
        shape=(len(barcodes), len(remo_ids)),
        dtype=np.int32
    )
    
    # Build var from metadata
    meta = metadata().to_pandas().set_index('REMO')
    var = meta.reindex(remo_ids)
    
    # Build obs
    obs = pd.DataFrame(index=pd.Index(barcodes, name='barcode'))
    
    logger.info(f'Created AnnData: {len(barcodes)} cells × {len(remo_ids)} modules')
    
    return AnnData(X=X, obs=obs, var=var)
