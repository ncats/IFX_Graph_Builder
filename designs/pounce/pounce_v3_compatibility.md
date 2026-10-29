# POUNCE v2/v3 workbook compatibility

## Source inspected

The v3 source set is the nine workbooks in
`~/Downloads/drive-download-20261027T191939Z-1-001.zip` (received 2026-10-27).
The v2 comparison set is `input_files/manual/pounce/OLD/templates/`.
The earlier `~/Downloads/updated again templates/` set has the same sheet and
header structure; one StatsResults definition changed in the final download.

## Format boundary

Detect each workbook by its role-specific sheet and map-column header. V2 maps
use `Submitter_Variable_Name`; v3 maps use `Submitter_Column_Name`. Filenames
are not trusted. Mixed headers within one workbook are unsupported.

`workbook_format.py` normalizes in memory before the shared metadata parser and
graph builder consume a workbook. The uploaded source file is preserved. The
normalizer records the original sheet names for validation reports.

| V3 input | Shared parser representation | Reason |
| --- | --- | --- |
| `BiospecimenMap`, `BiospecimenMeta` in proteomics Project | `BioSampleMap`, `BioSampleMeta` | Same project biosample and biospecimen row structure |
| `PeakDataMeta`, `PeakData` in proteomics Experiment | `ProteinDataMeta`, `ProteinData` | ProteinMap identifies the analyte family |
| `runbiosample_id` | `run_biosample_id` | Same run identity |
| `ncats_lab_groups` | `lab_groups` | Same NCATS lab-group concept |
| `exposureN_duration` | `exposureN_time` | Same duration field |
| `metabfeature_id` | `metab_id` for the measured row and matrix join | Feature identity, not a canonical metabolite ID |
| `standardized_metab_id` | separate parsed field and resolver input | Submitted standardized ID; may be blank |
| `ES_1`, `ES_1_Pval`, `ES_1_adjPval` | legacy primary effect-size description fields | Same primary comparison descriptions |

For v3 metabolites, the measured feature ID creates the `MeasuredMetabolite`
identity and matches PeakData and EffectSize rows. It is **not** sent to the
canonical metabolite resolver. Only submitted `standardized_metab_id` and
explicit alternate IDs are eligible for that edge. V2 keeps its existing
primary-ID behavior. No adapter-side cross-identifier reconciliation is added.

Common validation rules remain shared. The StatsResults required effect-size
fields are checked per workbook: v2 requires `ES`, `ESPval`, `ESadjPval`; v3
requires the primary `ES_1` description and `ES_1_prefix`. The v3
`ES_resultN`/`Pval_resultN`/`AdjPval_resultN` map keys are recognized and their
configured column names are checked against EffectSize. Their map keys define
the statistic type even when submitters use custom column prefixes.

Compatible mixed-version submissions are accepted. The project-to-experiment
biosample references, experiment-to-stats analyte mapping, and data matrix IDs
are validated from the normalized structures. A mismatch remains a validation
error tied to its submitted workbook. Validation messages translate normalized
sheet and key labels back to those in the uploaded file.

## Validation and known source issues

The existing POUNCE unit tests and new format tests cover detection, conflicting
headers, normalization, and v3 metabolite resolver semantics. The final v3
templates were parsed directly for metabolomics, proteomics, and transcriptomics;
the v2 transcriptomics template and a v2 Project plus v3 Experiment/StatsResults
set also parsed successfully.

The supplied v3 proteomics sample has internal reference errors: its
`BiospecimenMap` names an `exposures` column absent from `BiospecimenMeta`, and
its RunBioSampleMeta biosample IDs do not match Project biosample IDs. These
remain validation errors. Several supplied templates also leave required
content fields blank; support for the format does not treat those values as
valid submissions.

No Snakemake run or ETL database build was performed. The user should run the
targeted POUNCE ETL against a corrected v3 workbook set before a deployment.

## Submission policy update (2026-10-29)

Haley confirmed that each experiment may omit its StatsResults workbook. The
workbook and JSON interfaces keep an aligned slot for each experiment; a
missing StatsResults file is `null` in JSON. The parser and graph builder skip
that slot without shifting later stats files to a different experiment.

V3 ExperimentMeta omits `experiment_type`, so that field is optional there while
`platform_type` remains required. V2 still requires `experiment_type` and
checks its accepted values. The BioSampleMap
`disease_names` column mapping is optional; `N/A` row values were already
treated as missing data. The TMB lab-group label is accepted. `3D tissue` is
accepted as a biospecimen type because the field is stored as free text; the
supplied Excel template's accepted-values description should also be updated
by its owner.

GeneMeta omission and the reported MetabMap error remain open pending a data
contract decision and the failing workbook/error, respectively. The proposed
cross-category requirement for disease, exposure, condition, phenotype, or
demographic data is deferred.
