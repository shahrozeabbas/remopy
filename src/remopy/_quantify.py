'''Fragment quantification into REMO.'''

import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse
from anndata import AnnData

from remopy._data import modules, metadata


def quantify(
    fragments: str,
    min_fragments: int = 1000,
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
    keep = cell_counts.filter(pl.col('len') >= min_fragments).select('barcode')
    n_cells_before = cell_counts.height
    n_cells_after = keep.height
    frags = frags.join(keep, on='barcode')
    print(f'Filtered to {n_cells_after}/{n_cells_before} cells with >= {min_fragments} fragments')
    
    # Find overlaps between fragments and modules
    print('Finding overlaps...')
    overlaps = pb.overlap(frags, mods, how='inner', suffix='_mod')
    del frags
    
    # Count fragments per cell per module
    print('Counting fragments per module...')
    counts = overlaps.group_by(['barcode', 'REMO']).len().rename({'len': 'n'})
    del overlaps
    
    # Build sparse matrix
    print('Building count matrix...')
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
    
    print(f'Created AnnData: {len(barcodes)} cells × {len(remo_ids)} modules')
    
    return AnnData(X=X, obs=obs, var=var)
