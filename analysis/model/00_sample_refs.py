from utils_v10 import *

refs_path = sys.argv[1]
outdir = sys.argv[2]
sampling_size = int(sys.argv[3])  # 80
subsampling_n = int(sys.argv[4])  # 25
rng = np.random.default_rng(42)

refs = SeqIO.to_dict(SeqIO.parse(refs_path, "fasta"))
os.makedirs(f"results/{outdir}/sampled_contigs", exist_ok=True)
output_file = f"results/{outdir}/sampled_contigs/sampled_contigs_{subsampling_n}.fasta"

all_sampled_contigs = []
selected_sampling_contigs = {}
contigs_lens = []

for seq in refs.values():
    contigs, contigs_no = sampling(seq, sampling_size, rng)
    all_sampled_contigs.append(contigs)
    contigs_lens.append(contigs_no)

subsampled = subsample_contigs_dicts(all_sampled_contigs, subsampling_n, rng)

with open(output_file, "w") as f:
    SeqIO.write(subsampled.values(), f, "fasta")
