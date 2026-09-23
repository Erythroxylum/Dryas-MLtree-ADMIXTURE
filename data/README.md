# Input data

Place these files here before running the workflow:

- `s384_raw.vcf.gz`: full 384-sample ipyrad VCF with `GT:DP:CATG` fields.
- `Dryas_sampledata.csv`: metadata containing `sampleID`, `s380`, `Plate`, `geoID`, and `spID`.

The `s380` column identifies the 380 *Dryas* ingroup samples; the remaining four samples are retained as outgroups for the ML tree. Raw and cleaned VCFs are intentionally excluded from GitHub.

