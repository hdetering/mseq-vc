module load miniconda
source activate simtools

# global parameters
nclones=6
nsamples=5
nrep=10
outdir="sim_prep"
art="$HOME/src/art_bin_MountRainier/art_illumina"

for tt in us ms hs; do 
  for cvg in 50 100 300; do 
    python src/simtools/simtools.py scenario \
      --nclones $nclones \
      --nsamples $nsamples \
      --ttype $tt \
      --nrep $nrep \
      --out $outdir \
      --seq-read-gen true \
      --seq-coverage $cvg \
      --seq-art-path $art 2>&1
  done 
done > sim_prep.log
