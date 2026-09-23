# Dryas ML tree and ADMIXTURE

Reproducible filtering and reanalysis of the Pan-*Dryas* GBS SNP dataset after detection of a Plate 1/AK2019 genotype-calling artifact.

## Inputs

Start with the full, unfiltered 384-sample VCF exported by ipyrad (`s384_raw.vcf.gz`) and the sample metadata CSV. The VCF must contain `GT`, `DP`, and ipyrad `CATG` genotype fields; the metadata must contain `sampleID`, `s380`, `Plate`, `geoID`, and `spID`. The raw VCF is not stored here because it exceeds GitHub's file-size limit.

## Workflow

Create the software environment:

```bash
conda env create -f environment.yml
conda activate dryas-ml-admixture
```

Detect and remove the batch-associated loci:

```bash
python scripts/01_scan_vcf.py \
  --vcf data/s384_raw.vcf.gz \
  --metadata data/Dryas_sampledata.csv \
  --outdir results/batch_qc

python scripts/02_find_batch_loci.py --qc-dir results/batch_qc
python scripts/03_select_artifact_loci.py --qc-dir results/batch_qc

python scripts/04_filter_vcf.py \
  --vcf data/s384_raw.vcf.gz \
  --qc-dir results/batch_qc \
  --outdir results/cleaned_vcfs
```

The primary filter masks genotypes outside DP 15–150, masks heterozygotes with alternate-read balance outside 0.25–0.75, removes 637 batch-associated ipyrad loci, and recalculates MAC and call rate. It produced 95,479 SNPs in the 384-sample dataset and 24,160 SNPs in the 380-sample, ≥80%-called dataset.

Run two independently seeded ADMIXTURE replicates for K = 2–19 and plot the lower-CV replicate at each K:

```bash
bash scripts/05_run_admixture.sh \
  results/cleaned_vcfs/s380_batchclean_ab25_locusFDR05_mac4_dp15-150_miss80.vcf.gz \
  results/admixture

Rscript scripts/06_plot_admixture.R \
  results/admixture \
  data/Dryas_sampledata.csv \
  figures/admixture_cleaned.pdf \
  12 19
```

The LD-pruning window is explicitly 50 kb (`--indep-pairwise 50kb 5 0.2`). ADMIXTURE component colors are arbitrary and should not be interpreted as homologous among K values without component matching.

Build and plot the maximum-likelihood tree:

```bash
python scripts/07_vcf_to_iupac.py \
  --vcf results/cleaned_vcfs/s384_batchclean_ab25_locusFDR05_mac4_dp15-150_miss10.vcf.gz \
  --output results/ml_tree/dryas_cleaned.phy

bash scripts/08_run_iqtree.sh \
  results/ml_tree/dryas_cleaned.phy \
  results/ml_tree/dryas_cleaned

Rscript scripts/09_plot_tree.R \
  results/ml_tree/dryas_cleaned.treefile \
  data/Dryas_sampledata.csv \
  figures/dryas_cleaned_ml_tree.pdf
```

IQ-TREE uses `GTR+ASC` for the SNP-only alignment with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates.

## Figures

Final figures are written to `figures/` by the plotting scripts.

## References

Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research* 19:1655–1664.

Minh, B.Q. et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37:1530–1534.
