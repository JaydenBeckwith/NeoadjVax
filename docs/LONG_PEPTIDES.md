# Long peptides and the vaccine construct

pVACseq ranks 8–11mers for class I and 12–18mers for class II. A synthetic
long peptide vaccine is not made of those. It is made of 20–30mers that
carry the mutation with its natural flanking sequence, so the patient's
own antigen-presenting cells process and present it on whichever allele
fits — rather than the construct betting everything on one predicted
epitope and one predicted allele.

Turning the candidate table into that peptide list was the step still
being done by hand. `DESIGN_LONG_PEPTIDES` does it, and writes
`<patient>/long_peptides/<patient>.<class>.long_peptides.tsv` plus a
matching FASTA.

## Where the sequence comes from

No candidate table carries flanking sequence, so `GENERATE_PROTEIN_FASTA`
regenerates the mutant protein context with
`pvacseq generate_protein_fasta`, from the same VCF and the same phased
proximal-variant VCF pVACseq itself used. Proximal variants therefore
appear in the long peptide exactly as they appeared in the epitope.

## How the window is placed

1. **Find the mutation.** The MT and WT sequences of each fasta entry are
   compared residue by residue. This deliberately does not parse the fasta
   headers, whose field layout moves between pVACtools versions.
2. **Place the window.** A substitution is centred. A frameshift or
   insertion is *anchored at the junction*, because the novel stretch runs
   to the end of the protein and the part nearest the junction is what the
   immune system has not seen.
3. **Pull it inside the protein** at either terminus, and shrink it when
   the protein is shorter than `--long_peptide_max_length`.
4. **Trim awkward ends.** Residues in `--long_peptide_bad_nterm` (default
   `QNCM`) and `--long_peptide_bad_cterm` (default `CMP`) are trimmed off
   while the mutation keeps `--long_peptide_min_flank` residues of context
   on both sides and the peptide stays at least
   `--long_peptide_min_length`. Trimming never cuts past the mutation.
5. **Merge neighbours.** Two mutations close enough to fall inside one
   window become one peptide, both named in the `Covers` column.

`Contains Best Peptide` records whether the predicted epitope ended up
fully inside the designed peptide. It normally does; when it does not, the
peptide also carries a warning.

## Synthesis notes

The `Warnings` column is advisory and nothing is dropped for it:

- cysteine count (free cysteines oxidise; two or more can scramble)
- methionine count (oxidation)
- GRAVY hydropathy above `--long_peptide_max_gravy`, which tends to mean
  poor solubility
- a single-residue run at or above `--long_peptide_max_run`
- ambiguous or stop residues in the sequence

`Cys`, `Met`, `GRAVY`, `Net Charge` and `Longest Run` are reported as
columns too, so whoever places the order can apply their own synthesiser's
rules rather than these.

## The construct

`--vaccine_vector true` runs pVACvector on the designed peptide FASTA. It
orders the peptides into one string-of-beads construct and chooses spacers
so the junctions between peptides do not themselves create epitopes, which
is what you want for a DNA or RNA vaccine or a single long insert.

It is off by default because it searches over peptide orderings, rerunning
the binding predictors on every junction of every candidate ordering, and
that is a long job for something most runs do not need.

## Parameters

| parameter | default | what it does |
|---|---|---|
| `--long_peptides` | `true` | design peptides at all |
| `--long_peptide_max_length` | 25 | longest peptide to design |
| `--long_peptide_min_length` | 15 | shortest acceptable after trimming |
| `--long_peptide_min_flank` | 3 | context kept either side of the mutation |
| `--long_peptide_bad_nterm` | `QNCM` | trimmed from the N terminus where the flanks allow |
| `--long_peptide_bad_cterm` | `CMP` | trimmed from the C terminus where the flanks allow |
| `--long_peptide_top` | 20 | peptides per MHC class |
| `--vaccine_vector` | `false` | run pVACvector on the result |

Rejected candidates are skipped. Everything else is designed in
`Priority Rank` order, so `--long_peptide_top` takes the best ones.
