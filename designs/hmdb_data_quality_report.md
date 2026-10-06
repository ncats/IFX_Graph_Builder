# Draft HMDB data-quality report

Status: working draft for eventual submission to the HMDB maintainers.

This report collects possible source-data issues observed while integrating
HMDB into the NCATS metabolite harmonization workflow. It compares two
official HMDB 5.0 downloads and records discrepancies without assuming which
artifact should be corrected. Every proposed correction should be rechecked
against the current HMDB website before submission.

## Datasets examined

- `hmdb:metabolites_xml:5.0`: `hmdb_metabolites.zip`
- `hmdb:structures_sdf:5.0`: `structures.zip`
- HMDB version: 5.0
- Version date: 2021-11-17
- Download date used by NCATS: 2026-06-30

The XML contains 217,920 unique primary accessions and the SDF contains
217,776. Every SDF accession occurs in the XML; 144 accessions occur only in
the XML. Neither artifact contains duplicate primary accessions.

## 1. Chemistry-bearing records omitted from the SDF

Of the 144 XML-only accessions:

- 123 have formula, average mass, monoisotopic mass, and InChIKey;
- 120 also have SMILES;
- 119 also have InChI; and
- 21 have no chemistry fields in either examined artifact.

The XML-only chemistry appears internally reliable. All 123 formulas agree
with their reported monoisotopic masses within 0.00001 Da. All 120 SMILES
values parse with RDKit and reproduce the XML-reported InChIKey exactly.

Twenty-three of these omitted records are specific PI/PIP molecular species
mapped by HMDB to the broad KEGG record `C00626`:

| HMDB ID | HMDB name |
|---|---|
| HMDB0009834 | PI(18:1(9Z)/16:0) |
| HMDB0009837 | PI(18:1(9Z)/18:1(9Z)) |
| HMDB0009843 | PI(18:1(9Z)/20:3(8Z,11Z,14Z)) |
| HMDB0009871 | PI(20:1(11Z)/18:1(9Z)) |
| HMDB0009873 | PI(20:1(11Z)/20:4(5Z,8Z,11Z,14Z)) |
| HMDB0009874 | PI(20:1(11Z)/20:4(8Z,11Z,14Z,17Z)) |
| HMDB0009875 | PI(20:2(11Z,14Z)/16:0) |
| HMDB0009888 | PI(20:3(8Z,11Z,14Z)/18:1(11Z)) |
| HMDB0009891 | PI(20:3(8Z,11Z,14Z)/20:3(5Z,8Z,11Z)) |
| HMDB0009892 | PI(20:3(8Z,11Z,14Z)/20:3(8Z,11Z,14Z)) |
| HMDB0009893 | PI(20:4(5Z,8Z,11Z,14Z)/16:0) |
| HMDB0009895 | PI(20:4(5Z,8Z,11Z,14Z)/18:1(11Z)) |
| HMDB0009898 | PI(20:4(5Z,8Z,11Z,14Z)/20:1(11Z)) |
| HMDB0009900 | PI(20:4(8Z,11Z,14Z,17Z)/18:0) |
| HMDB0009904 | PI(20:4(8Z,11Z,14Z,17Z)/20:1(11Z)) |
| HMDB0009905 | PI(22:2(13Z,16Z)/16:0) |
| HMDB0009906 | PI(22:2(13Z,16Z)/18:2(9Z,12Z)) |
| HMDB0009909 | PI(22:3(10Z,13Z,16Z)/18:2(9Z,12Z)) |
| HMDB0009911 | PI(22:3(10Z,13Z,16Z)/18:3(9Z,12Z,15Z)) |
| HMDB0009916 | PI(22:5(4Z,7Z,10Z,13Z,16Z)/16:0) |
| HMDB0009918 | PI(22:5(7Z,10Z,13Z,16Z,19Z)/16:0) |
| HMDB0009927 | PIP(16:0/20:1(11Z)) |
| HMDB0009931 | PIP(16:0/20:4(5Z,8Z,11Z,14Z)) |

For example, XML record `HMDB0009911` supplies formula `C49H83O13P`, average
mass `911.1493`, monoisotopic mass `910.557129254`, InChIKey
`OMBBFDFDCMFWTB-YMBWUVAOSA-N`, SMILES, and InChI, but the accession is absent
from `structures.sdf`.

A machine-readable attachment containing all 144 XML-only accessions and the
123 chemistry-bearing records should be generated before submission.

## 2. Conflicting structures for the same accession

Thirty-nine shared accessions disagree on formula, InChI, InChIKey, and both
mass fields. Thirteen have different first InChIKey blocks, indicating
different molecular connectivity rather than a formatting or protonation
difference.

| HMDB ID | HMDB name | SDF formula | XML formula | SDF monoisotopic mass | XML monoisotopic mass |
|---|---|---:|---:|---:|---:|
| HMDB0013132 | Hydroxyvalerylcarnitine | C12H23NO5 | C12H25NO5 | 261.157622845 | 263.173272915 |
| HMDB0034051 | 6-Galloylglucose | C17H19BrN4O2 | C13H16O10 | 390.069138519 | 332.074346732 |
| HMDB0241915 | N-Nervonoyl Valine | C11H11NO10 | C29H55NO3 | 317.040489878 | 465.418194635 |
| HMDB0241916 | N-Oleoyl Glycine | C20H21FN6O5 | C20H37NO3 | 444.155746017 | 339.277344055 |
| HMDB0241917 | N-Linoleoyl Glycine | C17H33O2 | C20H35NO3 | 269.248055300 | 337.261693991 |
| HMDB0241918 | N-Arachidonoyl Glycine | C15H19NO9 | C22H35NO3 | 357.105981196 | 361.261693991 |
| HMDB0241919 | N-Palmitoyl Alanine | C21H34O5S | C19H37NO3 | 398.212695369 | 327.277344055 |
| HMDB0241920 | N-Palmitoyl Arginine | C21H34O6S | C22H44N4O3 | 414.207609989 | 412.341341293 |
| HMDB0241921 | N-Palmitoyl Asparagine | C10H12O6S | C20H38N2O4 | 260.035459280 | 370.283157712 |
| HMDB0241922 | N-Palmitoyl Aspartic acid | C26H43NO10S2 | C20H37NO5 | 593.232838937 | 371.267173295 |
| HMDB0241923 | N-Palmitoyl Cysteine | C14H22O4 | C19H37NO3S | 254.151809188 | 359.249415229 |
| HMDB0241924 | N-Palmitoyl Glutamine | C18H33NO4 | C21H40N2O4 | 327.240958547 | 384.298807776 |
| HMDB0242435 | Deoxycholylmethionine | C29H49NO5S | C15H20O3 | 523.333144856 | 248.141244504 |

The corresponding formulas, SMILES, InChIs, and InChIKeys should be included
in the submission attachment. Several SDF formulas appear incompatible with
the displayed HMDB name, but this report does not select a preferred artifact
without independent verification.

## 3. Protonation-state inconsistencies between XML and SDF

The other 26 structural-property disagreements preserve the first two
InChIKey blocks. In each case the XML formula generally has one additional
hydrogen, the XML InChI has a protonation layer, and the final InChIKey block
changes from `N` in the SDF to `O` in the XML. Most are acylcarnitines.

| HMDB ID | Name | SDF formula | XML formula |
|---|---|---:|---:|
| HMDB0000062 | L-Carnitine | C7H15NO3 | C7H16NO3 |
| HMDB0000201 | L-Acetylcarnitine | C9H17NO4 | C9H18NO4 |
| HMDB0000222 | Palmitoylcarnitine | C23H45NO4 | C23H46NO4 |
| HMDB0000736 | Isobutyryl-L-carnitine | C11H21NO4 | C11H22NO4 |
| HMDB0000756 | Hexanoylcarnitine | C13H25NO4 | C13H26NO4 |
| HMDB0000791 | Octanoylcarnitine | C15H29NO4 | C15H30NO4 |
| HMDB0000824 | Propionylcarnitine | C10H19NO4 | C10H20NO4 |
| HMDB0000848 | Stearoylcarnitine | C25H49NO4 | C25H50NO4 |
| HMDB0002250 | Dodecanoylcarnitine | C19H37NO4 | C19H38NO4 |
| HMDB0005065 | Oleoylcarnitine | C25H47NO4 | C25H48NO4 |
| HMDB0005066 | Tetradecanoylcarnitine | C21H41NO4 | C21H42NO4 |
| HMDB0006455 | Arachidonoylcarnitine | C27H45NO4 | C27H46NO4 |
| HMDB0006461 | Linoelaidylcarnitine | C25H45NO4 | C25H46NO4 |
| HMDB0006469 | Linoleyl carnitine | C25H45NO4 | C25H46NO4 |
| HMDB0006509 | Nervonyl carnitine | C31H59NO4 | C31H60NO4 |
| HMDB0013128 | Valerylcarnitine | C12H23NO4 | C12H24NO4 |
| HMDB0013130 | Glutarylcarnitine | C12H21NO6 | C12H22NO6 |
| HMDB0013133 | Methylmalonylcarnitine | C11H19NO6 | C11H20NO6 |
| HMDB0013329 | trans-2-Tetradecenoylcarnitine | C21H39NO4 | C21H40NO4 |
| HMDB0013334 | 9,12-Hexadecadienoylcarnitine | C23H41NO4 | C23H42NO4 |
| HMDB0061677 | O-Adipoylcarnitine | C13H23NO6 | C13H24NO6 |
| HMDB0240585 | cis-4-Decenoylcarnitine | C17H31NO4 | C17H32NO4 |
| HMDB0240588 | Myristoleoylcarnitine | C21H39NO4 | C21H40NO4 |
| HMDB0240664 | Benzoylcarnitine | C14H19NO4 | C14H20NO4 |
| HMDB0240665 | Lignoceroylcarnitine | C31H61NO4 | C31H62NO4 |
| HMDB0240666 | 3-Methyladipoylcarnitine | C14H25NO6 | C14H26NO6 |

These may both be chemically intentional representations, but the downloadable
artifacts should ideally document or consistently apply the chosen protonation
convention.

## 4. XML records without chemistry

Twenty-one primary accessions occur only in the XML and have no formula, mass,
SMILES, InChI, or InChIKey in the examined downloads:

`HMDB0189583`, `HMDB0234122`, `HMDB0240149`, `HMDB0240204`,
`HMDB0241963`–`HMDB0241974`, `HMDB0242236`, `HMDB0242317`,
`HMDB0242419`, `HMDB0242542`, and `HMDB0242613`.

This may be expected for provisional or incompletely curated records. It is
included so the maintainers can distinguish intentional incompleteness from
the 123 XML-only records that do contain complete chemistry.

## 5. Specific molecular species mapped to broad KEGG identifiers

HMDB maps 187 PI/PIP records directly to KEGG `C00626`, which represents a
broader phosphatidylinositol concept. The HMDB records include many distinct
formulas and structures. This may be a legitimate broad-class cross-reference,
but consumers should not interpret it as exact chemical identity.

The same pattern may occur for other lipid-class KEGG mappings. A complete
machine-readable inventory should be generated from the post-rebuild generic-
versus-specific validation rather than relying on this example alone.

## Overall consistency

The two HMDB artifacts agree very well outside the exceptions above:

- formula, InChI, and InChIKey agree for 217,737 of 217,776 shared records
  (99.9821%);
- both reported mass fields agree exactly for 217,717 of 217,756 numeric pairs
  (99.9821%); and
- 217,014 SMILES strings are textually identical, while another 723 differing
  strings retain the same InChIKey.

## Before submission

- Recheck the 13 connectivity conflicts against current HMDB record pages.
- Generate CSV attachments for all 144 XML-only records, the 39 cross-artifact
  conflicts, and broad KEGG mappings.
- Separate confirmed errors from intentional protonation conventions and
  broad-class cross-references.
- Include checksums or release metadata for both HMDB 5.0 archives.
- Add NCATS contact information and the IFX/ODIN software version used for the
  analysis.
