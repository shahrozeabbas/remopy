'''Fragment quantification into REMO modules.'''

import numpy as np
import polars as pl
from scipy import sparse
from anndata import AnnData

from remopy._data import modules, metadata


def quantify(
    fragments: str,
    min_fragments: int = 1000,
) -> AnnData:
    '''
    Quantify fragments into REMO modules.
    
    Counts the number of fragments overlapping each REMO module per cell.
    Each fragment is counted once (the read support column is ignored).
    
    Parameters
    ----------
    fragments
        Path to fragments file (10x format: chrom, start, end, barcode, [count]).
        Can be gzipped.
    min_fragments
        Minimum total fragments per cell to include in output.
        
    Returns
    -------
    AnnData
        AnnData object with cells as obs and REMO modules as var.
        X contains fragment counts (sparse integer matrix).
        var contains module metadata (CREs, Bases, GC_mean, CL, etc.).
        
    Examples
    --------
    >>> import remopy as remo
    >>> adata = remo.quantify('fragments.tsv.gz')
    >>> adata
    AnnData object with n_obs × n_vars = 5000 × 340069
    '''
    try:
        import polars_bio as pb
    except ImportError:
        raise ImportError(
            'polars-bio is required for quantification. '
            'Install with: pip install polars-bio'
        )
    
    # Load REMO modules
    mods = modules()
    
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
    keep = cell_counts.filter(pl.col('len') >= min_fragments)['barcode']
    frags = frags.filter(pl.col('barcode').is_in(keep))
    
    n_cells_before = cell_counts.height
    n_cells_after = keep.len()
    print(f'Filtered to {n_cells_after}/{n_cells_before} cells with >= {min_fragments} fragments')
    
    # Find overlaps between fragments and modules
    print('Finding overlaps...')
    overlaps = pb.overlap(frags, mods, how='inner', suffix='_mod')
    
    # Count fragments per cell per module
    print('Counting fragments per module...')
    counts = (
        overlaps
        .group_by(['barcode', 'REMO'])
        .len()
        .rename({'len': 'n'})
    )
    
    # Build sparse matrix
    print('Building count matrix...')
    barcodes = sorted(counts['barcode'].unique().to_list())
    remo_ids = sorted(mods['REMO'].unique().to_list())
    
    bc_idx = {b: i for i, b in enumerate(barcodes)}
    mod_idx = {m: i for i, m in enumerate(remo_ids)}
    
    rows = [bc_idx[b] for b in counts['barcode'].to_list()]
    cols = [mod_idx[m] for m in counts['REMO'].to_list()]
    data = counts['n'].to_list()
    
    X = sparse.csr_matrix(
        (data, (rows, cols)),
        shape=(len(barcodes), len(remo_ids)),
        dtype=np.int32
    )
    
    # Build var from metadata
    meta = metadata().to_pandas().set_index('REMO')
    var = meta.reindex(remo_ids)
    
    # Build obs
    obs = pl.DataFrame({'barcode': barcodes}).to_pandas().set_index('barcode')
    
    print(f'Created AnnData: {len(barcodes)} cells × {len(remo_ids)} modules')
    
    return AnnData(X=X, obs=obs, var=var)
