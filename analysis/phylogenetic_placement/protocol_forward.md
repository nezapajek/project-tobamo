
# make a reference multiple alignment and make a tree
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --threads 32

# make a refernce multiple alignment, trimm and make a tree
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --trimal \
  --threads 32

# create query alignemnts and make phylogenetic placements
# untrimmed ref tree
python scripts/phylo_placement.py \
  --from-preprocess results/untrimmed/preprocess_summary.json \
  --contig-fasta data/tobamo_65_representative_contigs_wrong_orientation_corrected.fasta \
  --output-dir results/untrimmed/placement_65

## trimmed ref tree
python scripts/phylo_placement.py \
  --from-preprocess results/trimmed/preprocess_summary.json \
  --contig-fasta data/tobamo_65_representative_contigs_wrong_orientation_corrected.fasta \
  --output-dir results/trimmed/placement_65

# make a tree with each of the query sequences represented as a pendand edge
# untrimmed ref tree
gappa examine graft \
  --jplace-path results/untrimmed/placement_65/epa_output/epa_result.jplace \
  --out-dir results/untrimmed/gappa_65

# untrimmed ref tree
gappa examine graft \
  --jplace-path results/trimmed/placement_65/epa_output/epa_result.jplace \
  --out-dir results/trimmed/gappa_65