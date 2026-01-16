# remopy

Python implementation of [REMO.v1.GRCh38](https://github.com/stuart-lab/REMO.v1.GRCh38), the R data package from the [Stuart Lab](https://github.com/stuart-lab).

REMO (Regulatory Element MOdules) provides pre-defined, cell-type annotated regulatory element groupings for single-cell chromatin accessibility analysis.

## Installation

```bash
# Core data package (just polars)
pip install remopy

# With fragment quantification support
pip install remopy[quantify]
```

## Quick Start

### Data Access

```python
import remopy as remo

# Load module coordinates (1.5M CRE intervals → 340k modules)
modules = remo.modules()
print(modules.head())

# Load module metadata
metadata = remo.metadata()
print(metadata.columns)  # ['REMO', 'CREs', 'Bases', 'Chromosome', 'GC_mean', 'CL']

# Get modules associated with a cell type
terms = remo.terms()
t_cell_modules = terms.get('T cell', [])

# Get cell types present in a tissue
tissues = remo.tissues()
brain_cell_types = tissues.get('Brain', [])
```

### Fragment Quantification (scATAC-seq)

Skip peak calling entirely — quantify fragments into REMO:

```python
import scanpy as sc
import remopy as remo

# Quantify fragments into modules (requires polars-bio)
adata = remo.quantify('fragments.tsv.gz', min_fragments=1000)

# Standard scanpy workflow
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

## Data Contents

| Data | Description |
|------|-------------|
| `modules()` | 1,507,327 CRE intervals grouped into 340,069 modules |
| `metadata()` | Module-level stats: CRE count, bases, GC content, cell ontology |
| `terms()` | Cell type name → module ID mappings (144 cell types) |
| `ontology()` | Cell Ontology ID → module ID mappings |
| `tissues()` | Tissue → cell type mappings (25 tissues) |

## Why REMO?

- **No peak calling needed**: Use pre-defined, validated features
- **Reproducible**: Same features across all datasets
- **Cell-type annotated**: Modules linked to Cell Ontology terms
- **Fast**: Direct fragment → module quantification

## Citation

Lim C, et al. Regulatory element modules as universal features for single-cell chromatin analysis. (2025)

[Preprint on bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.10.692786v1)

## License

Artistic License 2.0
