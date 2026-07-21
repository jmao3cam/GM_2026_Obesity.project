# RRBS reference files

Reference/annotation files used by `RRBS/01_data_loading.R` (MspI + TSS filters).

## `hg38_TSS.bed`

This is actually a **gzip-compressed UCSC `refGene.txt` table** with a `.bed`
extension (confirmed from the file's magic bytes — it's not plaintext BED).
That's fine as-is: `data.table::fread()` (used in `01_data_loading.R`)
auto-detects gzip regardless of file extension, so no change needed. Kept
compressed here since it's ~8.6MB compressed and comfortably under GitHub's
size limits.

## `build_mspi_fragments.py`


- MspI cuts `CCGG` between the two Cs (`C^CGG`) — standard, well-documented
  enzyme behaviour, not inferred from your file.
- The specific cutting/boundary logic (fragments run edge-to-edge across
  each chromosome, first fragment starts at 0, last fragment ends at the
  chromosome's full length) was confirmed by checking that chrM's last
  fragment ends at exactly 16,569 — hg38's true chrM length — and that
  chr1's first fragment (0–10,482) matches hg38's known leading assembly
  gap.
- Only standard chromosomes (chr1–22, X, Y, M) appear in your file, so the
  script defaults to that set.


### generating `MspI_fragments_hg38.bed`

```bash
pip install pyfaidx --break-system-packages
python build_mspi_fragments.py --fasta /path/to/hg38.fa --out MspI_fragments_hg38.bed
```
