# Draft RefMet data-quality report

Status: working draft for eventual submission to the RefMet maintainers.

This report collects possible source-data issues observed while integrating
RefMet into the NCATS metabolite harmonization workflow. It should be reviewed
and supplemented with the post-rebuild validation results before submission.

## Dataset examined

- Registry artifact: `refmet:metabolites_csv:sha256-ef152c85ded4`
- RefMet download/version date: 2026-09-09
- Rows: 208,170
- Fields used here: `refmet_id`, `refmet_name`, `formula`, `exactmass`,
  `inchi_key`, and external database identifiers

## 1. Formula and exact-mass discrepancies

We compared each populated RefMet `exactmass` with the neutral monoisotopic
mass calculated from its `formula`. Calculations used RDKit periodic-table
most-common-isotope masses and included all dot-separated formula components.

The fields are highly consistent overall:

- 208,141 records had both values and could be compared;
- all 208,141 populated formulas were parseable;
- median absolute difference was 0.000001526 Da;
- 99th-percentile absolute difference was 0.000025326 Da; and
- 18 records differed by more than 0.5 Da.

The following are the 18 exceptions. `Delta` is reported exact mass minus the
mass calculated from the formula.

| RefMet ID | Name | Formula | Reported exact mass | Calculated mass | Delta (Da) | Reference ID |
|---|---|---:|---:|---:|---:|---|
| RM0004956 | Methylbromoeudistomin D | C12H8Br2N2O | 355.898291 | 353.900337 | +1.997954 | — |
| RM0005392 | Tetrachlorocyclohexene | C6H6Cl4 | 219.919411 | 217.922361 | +1.997050 | — |
| RM0021730 | 3,4,5,6-Tetrachlorocyclohexene | C6H6Cl4 | 219.919411 | 217.922361 | +1.997050 | PubChem 15691; ChEBI 28988; KEGG C06989 |
| RM0108748 | 9-Methyl-7-bromoeudistomin D | C12H8Br2N2O | 355.898291 | 353.900337 | +1.997954 | PubChem 129958; ChEBI 34514; KEGG C13748 |
| RM0108932 | Oxiconazole nitrate | C18H14Cl4N4O4 | 491.973966 | 489.976916 | +1.997050 | PubChem 9556529; ChEBI 7826; KEGG C08075 |
| RM0118236 | Mitobronitol | C6H12Br2O4 | 307.908187 | 305.910233 | +1.997954 | PubChem 656655; ChEBI 34853; KEGG C13522 |
| RM0137167 | Miconazole nitrate | C18H15Cl4N3O4 | 478.978717 | 476.981667 | +1.997050 | PubChem 68553; ChEBI 82892; KEGG C08070 |
| RM0137333 | IAARh123 | C28H19IN5O5 | 630.043249 | 632.043092 | -1.999843 | PubChem 656577; ChEBI 5845; KEGG C11599 |
| RM0157487 | Glcbeta1-6Glcbeta-Caldarchaeol | C98H192O16 | 1626.424395 | 1625.421040 | +1.003355 | PubChem 11973051; ChEBI 34773; LipidMaps LMGL02060005; KEGG C13874 |
| RM0159336 | IACI | C26H27IN6O3 | 596.119094 | 598.118937 | -1.999843 | PubChem 656579; ChEBI 5846; KEGG C11600 |
| RM0189207 | N-Acetylamylamine | C7H15NO | 113.120449 | 129.115364 | -15.994915 | PubChem 221432 |
| RM0189210 | N-Propionylcysteine | C6H11NO3S | 145.073894 | 177.045964 | -31.972070 | PubChem 151048 |
| RM0189229 | N-Propionylamylamine | C8H17NO | 127.136099 | 143.131014 | -15.994915 | PubChem 4088627 |
| RM0189231 | N-Propionylhistamine | C8H13N3O | 151.110947 | 167.105862 | -15.994915 | PubChem 24284480 |
| RM0189233 | N-Propionylputrescine | C7H16N2O | 128.131348 | 144.126263 | -15.994915 | PubChem 28754146 |
| RM0189234 | N-Propionyltaurine | C5H11NO4S | 149.068809 | 181.040879 | -31.972070 | PubChem 53723115 |
| RM0200310 | Triforin | C10H14Cl6N4O2 | 433.921842 | 431.924792 | +1.997050 | PubChem 33565; ChEBI 9715; KEGG C10960 |
| RM0204720 | 6-Bromo-3'-methylflavone | C16H11BrO2 | 314.994000 | 313.994242 | +0.999758 | PubChem 688751 |

Patterns worth checking at the source include four records whose reported mass
is lower by the mass of one oxygen atom, two lower by the mass of one sulfur
atom, and several halogen-containing compounds whose reported value is about
one or two daltons above the formula-derived monoisotopic mass. These may mix
formula errors, exact-mass errors, and isotope conventions; the table reports
the observations without choosing which field should be changed.

## 2. Malformed InChIKey

| RefMet ID | Name | Reported InChIKey | Expected from PubChem cross-reference | Reference |
|---|---|---|---|---|
| RM0139008 | 2-Dodecylbenzenesulfonic acid | `WBINGBSDOWNP-UHFFFAOYSA-N` | `WBIQQQGBSDOWNP-UHFFFAOYSA-N` | PubChem 25457; ChEBI 149774; HMDB0031031 |

The first block of the reported value has 12 characters rather than the 14
required by the InChIKey format. The expected value above is included as a
candidate correction and should be independently confirmed before submission.

### Duplicate InChIKeys with incompatible chemistry

The file contains 12 InChIKeys used by more than one RefMet record. Most are
duplicate names, stereochemical variants, or differ only in harmless decimal
precision. Two full InChIKeys, however, are assigned to records with different
formulas and masses:

| Reported InChIKey | RefMet ID | Name | Formula | Exact mass |
|---|---|---|---|---:|
| `GXPOZJZFGYCEMT-VEKMRQFKSA-N` | RM0183628 | LPG 20:4(5Z,8Z,11Z,14Z)/0:0 | C26H45O9P | 532.280123 |
| `GXPOZJZFGYCEMT-VEKMRQFKSA-N` | RM0200938 | PG 20:4(5Z,8Z,11Z,14Z)/0:0 | C26H43O10P | 546.259388 |
| `SULIDBRAXVDKBU-PTGWMXDISA-N` | RM0174010 | LPC 0:0/18:1(9Z) | C26H52NO7P | 521.348142 |
| `SULIDBRAXVDKBU-PTGWMXDISA-N` | RM0133521 | PC 0:0/18:1(9Z) | C26H50NO8P | 535.327400 |

Because the complete InChIKey is identical while the oxygen count and exact
mass differ, at least one key in each pair is inconsistent with the associated
record. The two other duplicate-key groups flagged by an exact string audit
have the same formula and differ only in the last few reported mass decimals.

Nineteen connectivity blocks occur with more than one formula. Several are
expected acid/base or protonation pairs, so connectivity-block disagreement is
only a review signal rather than a confirmed-error count. The two full-key
cases above are the high-confidence subset.

## 3. Suspected incorrect cross-references

These mappings connected records with strongly incompatible identities or
masses in the harmonized graph. They need a final source-by-source verification
before being presented as requested corrections.

| RefMet ID | Suspected mapping problem | Current evidence |
|---|---|---|
| RM0139026 | ChEBI 49415 conflicts with PubChem 118701016 | ChEBI record is cobalt(3+); PubChem record is vitamin B12 |
| RM0200998 | ChEBI 147138 conflicts with PubChem 33572 | Linked records have masses approximately 2591 and 287.4 Da |
| RM0030757 | ChEBI 138856 conflicts with PubChem 18972857 | ChEBI record is oxolinic acid; linked masses are approximately 261.2 and 114.1 Da |
| RM0228128 | ChEBI 152031 conflicts with PubChem 118658 | Linked records have masses approximately 2661 and 158.3 Da |
| RM0227083 | ChEBI 110006 conflicts with HMDB0006044 | Linked records have masses approximately 516.7 and 165.2 Da |
| RM0228119 | HMDB0032386 conflicts with PubChem 11339 | Linked records have masses approximately 134.2 and 422.6 Da |

### Systematic source-file cross-reference audit

The RefMet file's own formula, exact mass, and reported InChIKey were compared
with the referenced PubChem and ChEBI records currently ingested by NCATS. This
is a candidate-generation audit: salts, protonation states, stereochemistry,
and generic structures can legitimately differ and must be reviewed before
calling an individual mapping wrong.

| Comparison | RefMet xref rows | Formula differs | InChIKey connectivity differs | Mass differs by >0.01 Da |
|---|---:|---:|---:|---:|
| PubChem | 35,471 | 378 | 233 | 143 |
| ChEBI | 25,050 | 997 of 24,996 comparable | 371 | 926 of 24,857 comparable |

Many of the largest new graph warnings are clear ChEBI cross-reference errors,
not errors in the RefMet chemistry fields. Representative examples are:

| RefMet record | RefMet chemistry | Referenced ChEBI record | Assessment |
|---|---|---|---|
| RM0157120, NA-Ser 7:0 | C10H19NO4; 217.131409; `OQVZMSFZPLEVSG-QMMMGPOBSA-N` | CHEBI:186510, PI(20:1/22:6); C51H85O13P; 936.57278 | Unrelated identity; RefMet agrees with PubChem 56989140 |
| RM0200095, oxygenated polyol | C12H26O5; 250.178025; `WMYINDVYGQKYMI-UHFFFAOYSA-N` | CHEBI:185935, PI(15:1/22:6); C46H75O13P; 866.49453 | Unrelated identity; RefMet agrees with PubChem 90038 |
| RM0154558, NA-Leu 10:0 | C16H31NO3; 285.230394; `KDQCSJBLQLMYLH-AWEZNQCLSA-N` | CHEBI:185722, PI(P-20:0/22:0); C51H99O12P; 934.68742 | Unrelated identity; RefMet agrees with PubChem 13890321 |
| RM0189227, N-Propionyl-DOPA | C12H15NO5; 253.095024; `AQESIQPVORIUAV-QMMMGPOBSA-N` | CHEBI:186338, PI(P-16:0/17:2); C42H77O12P; 804.51526 | Unrelated identity; RefMet agrees with PubChem 61152483 |
| RM0199360, Diofenolan | C18H20O4; 300.136160; `ZDOOQPFIGYHZFV-UHFFFAOYSA-N` | CHEBI:184766, PI(P-16:0/22:6); C47H79O12P; 866.53091 | Unrelated identity; RefMet agrees with PubChem 4578030 |

The PubChem audit also reveals a likely shifted identifier series. RefMet's
`NA-Arg` records from RM0154333 onward repeatedly point to a PubChem record one
carbon shorter than the RefMet name, formula, mass, and InChIKey. For example,
RM0154333 `NA-Arg 4:0` is C10H20N4O3 at 244.153541 Da but points to PubChem
171116623, which is C9H18N4O3 at 230.13789045 Da; PubChem 171116624 matches the
next formula in the series. This pattern should be reviewed as a batch rather
than as isolated discrepancies.

## 4. Specific lipids mapped to broader KEGG records

Some RefMet lipid molecular species share a KEGG compound identifier even
though their formulas and masses differ substantially. This may reflect a
legitimate broad-class cross-reference rather than an erroneous RefMet value,
but such mappings should not be interpreted as chemical identity.

Examples to review include KEGG C00958, C05973, C04438, C15647, C03819,
C18126, C18125, and C02530. The post-rebuild results should be used to expand
this section into a complete machine-readable attachment rather than relying
on examples alone.

## 5. Post-rebuild final-stage MW warnings

The completed `RaMP-ish - after more` run
(`20260928T140512Z-08c4fe77`) had 41 molecular-weight spread warnings in its
final stage (`stage-07-2ed3c998ab397b22`), down from 130 before final cleanup.
Thirty-three warning cliques contained 34 RefMet records; eight warning
cliques contained no RefMet record.

We reconstructed the effective InChIKey matching used by the stage in addition
to inspecting retained source cross-reference edges. This matters because
InChIKey unions are not stored as synthetic evidence edges. The reconstruction
used the reported key when present, otherwise the structure-derived key, and
used the first key block at or above the 500-Da cutoff and the first two blocks
below it.

### RefMet records with internally conflicting identity fields

These eight records are the strongest source-row correction candidates. In
each case, multiple external records agree on one chemistry while RefMet's
formula/mass, reported InChIKey, or both describe another. The first three
lipid cases are not merely an artifact of the 500-Da connectivity-block rule:
at least one external record carries the same complete reported InChIKey.

| RefMet ID | Name | RefMet formula / exact mass | Conflicting evidence | Assessment |
|---|---|---|---|---|
| RM0233955 | PG 18:2(9Z,12Z) | C24H43O10P / 522.2594 | Reported full key and both xrefs identify C44H79O10P / 798.541 | Formula/mass/name conflict with key and xrefs |
| RM0233954 | PS 18:2(9Z,12Z) | C24H42NO10P / 535.2546 | Reported full key and HMDB, LipidMaps, and PubChem peers identify C44H78NO10P / 811.536 | Formula/mass/name conflict with key and xrefs; clique is joined by InChIKey |
| RM0233963 | PI 18:2(9Z,12Z) | C27H47O13P / 610.2754 | Reported key and all three xrefs identify C47H83O13P / 886.557 | Formula/mass/name conflict with key and xrefs |
| RM0200497 | Carene | C1H160 / 173.252 | Key and ChEBI, HMDB, and PubChem agree on C10H16 / 136.125 | Apparent digit-transposition typo in formula with correspondingly wrong mass |
| RM0224683 | 2-Carene | C1H160 / 173.252 | Key and ChEBI, HMDB, and PubChem agree on C10H16 / 136.125 | Apparent digit-transposition typo in formula with correspondingly wrong mass |
| RM0188995 | 1,3-Diisopropylbenzene | C1H182 / 195.42415 | Key and ChEBI, HMDB, and PubChem agree on C12H18 / 162.141 | Apparent digit-transposition typo in formula with correspondingly wrong mass |
| RM0138911 | 3beta,7alpha-Dihydroxy-5-cholestenoic acid | C27H44O4 / 432.32396 | HMDB:HMDB0012454 and KEGG:C17335 agree with the C27 formula/mass; the reported key plus LipidMaps:LMST04010217, PubChem:3082147, and ChEBI:196598 identify the shorter C24 compound (C24H38O4 / 390.277) | Retain the HMDB and KEGG mappings; correct the reported key and remove or replace the LipidMaps, PubChem, and ChEBI mappings |
| RM0006550 | D-Ribulose | C5H10O6 / 166.04774 | Key, name, and ChEBI/PubChem/HMDB xrefs identify C5H10O5 / 150.0528 | Formula has one extra oxygen and mass follows the incorrect formula |

For RM0138911, likely C27 replacements include ChEBI:81015 and
LipidMaps:LMST04030241. The latter specifies the 25R stereoisomer and links
PubChem:77461123, whereas the HMDB structure has a different stereochemical
InChIKey and links PubChem:3081084. RefMet's name does not specify that side-chain
stereocenter, so selecting a replacement full InChIKey and PubChem record needs
a stereochemistry policy decision; the current C24 connectivity block
`PXHCARRJGFGPAC` is unambiguously wrong.

### High-confidence incompatible ChEBI cross-references

For these records, RefMet's name, formula, mass, and any non-ChEBI chemistry
evidence agree, while the RefMet-supplied ChEBI identifier points to a
different lipid or compound. These mappings account for 11 final-stage warning
cliques.

| RefMet ID | RefMet record | Suspected incompatible ChEBI xref |
|---|---|---|
| RM0045552 | GlcCer 14:1;O2(4E)/22:1;O; 741.575 Da | CHEBI:185532, PI(22:1/20:0); 948.667 Da |
| RM0226732 | alpha-Tocopherol succinate; 530.397 Da | CHEBI:22470; C29H50O2 / 430.381 Da |
| RM0174443 | DG 18:3/20:3/0:0; 640.507 Da | CHEBI:89189; C35H62O5 / 562.460 Da |
| RM0175770 | TG 16:0/18:1/20:3; 882.768 Da | CHEBI:187594, PI(13:0/18:3); 790.463 Da |
| RM0031923 | PI 16:1/11:2;O; 750.396 Da | CHEBI:88563, PI(16:1/18:1); 834.526 Da |
| RM0174446 | DG 18:3/15:0/0:0; 576.475 Da | CHEBI:89162, DG(18:3/20:3/0:0); 640.507 Da |
| RM0173312 | PA O-20:0/22:4; 766.588 Da | CHEBI:184961, PI(P-20:0/16:0); 850.594 Da |
| RM0163893 | PI 20:3/11:2;O; 802.427 Da | CHEBI:89131, PI(20:3/18:1); 886.557 Da |
| RM0174231 | DG 14:1/18:2/0:0; 562.460 Da | CHEBI:88721, DG(14:1/22:1/0:0); 620.538 Da |
| RM0098147 | PE 38:0/4:1; 829.656 Da | CHEBI:189610, PE(37:5); 751.515 Da |
| RM0199658 | [8]-Shogaol; C17H24O3 / 276.173 Da | CHEBI:174643, [8]-Shogaol labeled C19H28O3 / 304.204 Da |

### Representation differences requiring review, not automatic correction

Seven warning cliques reflect free-base/acid versus salt, counterion, nitrate,
or multicomponent representations. RefMet's mass is generally internally
consistent in these cases, but its xrefs may mix representation levels:

- RM0129769 Vasopressin is 1055.432 Da and agrees with HMDB and PubChem;
  CHEBI:9937 reports a dot-separated doubled/mixed formula at 2138.870 Da.
- RM0224347 6-Deoxyfagomine is the 131.095-Da free base, while its PubChem xref
  is the 167.071-Da hydrochloride.
- RM0108900 Naphazoline is the 210.116-Da free base, while CHEBI:7470 is the
  246.092-Da hydrochloride.
- RM0137167 Miconazole nitrate and RM0200855 Miconazole are distinct nitrate
  and base records joined through shared source mappings; their own masses are
  appropriate to their formulas.
- RM0161349 Delphinidin mixes chloride-form formula/mass/key evidence with
  neutral-form ChEBI and HMDB xrefs.
- RM0043113 Nylidrin has an HCl-form formula/mass but a neutral-form InChIKey
  and neutral ChEBI, HMDB, and PubChem xrefs.
- RM0200635 Berberrubine is the 322.108-Da cation/free form, while CHEBI:175111
  contains chloride and is 357.077 Da.

### RefMet present but not responsible for the warning

Seven RefMet records agree with their formula-derived mass, InChIKey peers,
and direct primary xrefs. Another member or mapping elsewhere in each clique
causes the warning: RM0153612 (arachidonic acid), RM0061038
(plastoquinol-1), RM0135878 (dopamine), RM0159410 (terfenadine), RM0136449
(nabilone), RM0126063 (pregnanediol), and RM0161190
(N-desmethylcitalopram).

## Before submission

- Recheck every proposed correction against the current RefMet website.
- Add the full post-rebuild list of incompatible external cross-references.
- Separate confirmed errors from valid broad-class mappings.
- Attach a CSV containing the affected RefMet IDs, fields, observed values,
  proposed values, evidence URLs, and validation status.
- Add NCATS contact information and the IFX/ODIN software version used for the
  analysis.
