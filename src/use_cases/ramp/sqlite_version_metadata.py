"""Human-facing release labels; exact Registry identities remain in provenance."""
import re

SOURCES = {
    'hmdb': ('HMDB', 'https://hmdb.ca/'),
    'chebi': ('ChEBI', 'https://www.ebi.ac.uk/chebi/'),
    'uniprot': ('UniProt', 'https://www.uniprot.org/'),
    'rhea': ('Rhea', 'https://www.rhea-db.org/'),
    'reactome': ('Reactome', 'https://reactome.org/'),
    'wikipathways': ('WikiPathways', 'https://www.wikipathways.org/'),
    'lipidmaps': ('Lipid Maps', 'https://www.lipidmaps.org/'),
    'refmet': ('RefMet', 'https://www.metabolomicsworkbench.org/databases/refmet/browse.php'),
    'pubchem': ('PubChem', 'https://pubchem.ncbi.nlm.nih.gov/'),
    'expasy': ('ExPASy ENZYME', 'https://enzyme.expasy.org/'),
    'pfocr': ('Pathway Figure OCR', 'https://pfocr.wikipathways.org/'),
}
DATASETS = {
    'ontology_full': 'Ontology', 'three_star_sdf': 'Structures',
    'metabolites_xml': 'Metabolites', 'proteins_xml': 'Proteins',
    'structures_sdf': 'Structures', 'cid_molecular_info': 'Molecular properties',
}


def release_label(source, dataset):
    version = str(dataset.get('version') or '')
    if version.startswith(('sha256-', 'deps-')):
        kind = 'Derived dataset' if version.startswith('deps-') else 'Content snapshot'
        date = dataset.get('download_date')
        return f'{kind}; downloaded {date}' if date else f'{kind}; release not recorded'
    if not version:
        return 'Release not recorded'
    if re.fullmatch(r'\d{4}-\d{2}-\d{2}', version):
        label = version
    elif source in ('rhea', 'chebi'):
        label = f'Release {version}'
    elif source in ('hmdb', 'reactome'):
        label = version if version.startswith('v') else f'v{version}'
    else:
        label = version
    date = dataset.get('version_date')
    return f'{label} ({date})' if date and date != version else label


def display_metadata(source, datasets):
    name, url = SOURCES.get(source, (source, ''))
    labels = {release_label(source, d) for d in datasets}
    if len(labels) == 1:
        version = next(iter(labels))
    else:
        version = '; '.join(sorted({
            f"{DATASETS.get(d['dataset'], d['dataset'].replace('_', ' ').capitalize())}: {release_label(source, d)}"
            for d in datasets
        }))
    return dict(data_source_name=name, data_source_url=url, data_source_version=version)
