# Representatives generation workflows

This README contains two clearly separated workflows:
1. **Virgaviridae representatives selection** (target class references)
2. **Non-virga representatives** used as model controls (neither of the two classes)

---

## PART 1 — Virgaviridae representatives selection

Input directory:
`/home/tobamo/analize/project-tobamo/analysis/data/domain_sci_input/ncbi_virgaviridae_20250120`

Main steps:
- Split/process sequences by genus: furovirus, goravirus, hordeivirus, pecluvirus, pomovirus, tobamovirus, tobravirus.
- Remove duplicate sequences per genus.
- Compute pairwise sequence identity within each genus (MAFFT-based alignment).
- Convert similarity to distance (`1 - identity`).
- Perform hierarchical clustering (single linkage).
- Select candidate representatives (longest sequence per cluster at identity threshold).
- Save outputs to `analysis/references/results/virga`:
	- `pairwise_aln/`
	- `clusters/`
	- `selected_species/`
	- `dendrograms/`
	- `dendrograms_selected_species/`

### Environment

`tobamo-model` is defined in:
`/home/tobamo/analize/project-tobamo/analysis/model/analysis_environment.yml`

```bash
conda env create -f ../model/analysis_environment.yml
conda activate tobamo-model
```

Note: `mafft` must be available in `PATH`.

### Run
From `analysis/references`:
```bash
python get_genus_representatives.py
```

Optional example:
```bash
python get_genus_representatives.py --identity-threshold 0.90 --max-workers 36 --overwrite
```

Final curation note:
The final selection decision was made by the domain scientist, and the curated selection is in:
`/home/tobamo/analize/project-tobamo/analysis/data/tobamo/reference_database.xlsx`

---

## PART 2 — Non-virga representatives (model controls)

Purpose:
Select non-virga control representatives for the model (neither of the two classes) for testing both the effectiveness of the similarity-based prefilter and the model’s ability to reject non-Virgaviridae sequences retained due to conserved domains or misannotations. 

Input FASTA (1,398 NCBI viral representatives (accessed January 16, 2025) duplicated at the genus level):
`/home/tobamo/analize/project-tobamo/analysis/data/non-virga_representatives/ncbi_virus_refseq_30k_dedupl_seqs_20250116.fasta`

Filtering happens at the DIAMOND TPDB2 step, so run only this part via script.

### Environment
```bash
conda create -n diamond python=3.9 -y
conda activate diamond
conda install -c bioconda diamond=0.9.14 biopython=1.81 pandas=2.1.1 megan=6.21.7 -y
```

### Run
From `analysis/references`:
```bash
python run_diamond_pipeline.py
```