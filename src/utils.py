import pandas as pd
from functools import lru_cache

GENE_ANNOTATIONS = {
    "mouse": "M023",
    "human": "G026",
    "Mus musculus": "M023",
    "Homo sapiens": "G026",
}


@lru_cache(maxsize=10)
def download_sample_metadata(srp: str, species: str = "mouse") -> pd.DataFrame:
    """
    Download and parse sample metadata for a specific SRA project.

    Parameters
    ----------
    srp : str
        SRA project ID (e.g., 'SRP115307')
    species : str, optional
        Species name, either 'mouse' or 'human' (default: 'mouse')

    Returns
    -------
    pd.DataFrame
        DataFrame containing sample metadata with columns specific to the SRA project

    Raises
    ------
    AssertionError
        If species is not in the supported GENE_ANNOTATIONS
    """

    assert species in GENE_ANNOTATIONS
    url = (
        f"https://duffel.rail.bio/recount3/{species}/data_sources/sra/"
        f"metadata/{srp[-2:]}/{srp}/sra.sra.{srp}.MD.gz"
    )
    sample_anno = pd.read_csv(url, delimiter="\t")
    return sample_anno


@lru_cache(maxsize=10)
def download_project_metadata(srp: str, species: str = "mouse") -> pd.DataFrame:
    """
    Download and parse project-level metadata for a specific SRA project.

    Parameters
    ----------
    srp : str
        SRA project ID (e.g., 'SRP115307')
    species : str, optional
        Species name, either 'mouse' or 'human' (default: 'mouse')

    Returns
    -------
    pd.DataFrame
        DataFrame containing project metadata with columns specific to the SRA project

    Raises
    ------
    AssertionError
        If species is not in the supported GENE_ANNOTATIONS
    """
    assert species in GENE_ANNOTATIONS
    url = (
        f"https://duffel.rail.bio/recount3/{species}/data_sources/sra/metadata/"
        f"{srp[-2:]}/{srp}/sra.recount_project.{srp}.MD.gz"
    )
    sample_anno = pd.read_csv(url, delimiter="\t")
    return sample_anno


@lru_cache(maxsize=10)
def download_gene_metadata(species: str = "mouse") -> pd.DataFrame:
    """
    Download and parse gene annotations for a specific species.

    Parameters
    ----------
    species : str, optional
        Species name, either 'mouse' or 'human' (default: 'mouse')

    Returns
    -------
    pd.DataFrame
        DataFrame containing gene annotations with columns:
        - seqname: chromosome or scaffold name
        - source: annotation source
        - feature: feature type
        - start: start position
        - end: end position
        - score: score (usually '.')
        - strand: strand (+/-)
        - frame: reading frame
        - gene_id: unique gene identifier
        - gene_type: type of gene
        - gene_name: common gene name

    Raises
    ------
    KeyError
        If species is not in the supported GENE_ANNOTATIONS
    """
    try:
        version = GENE_ANNOTATIONS[species]
    except KeyError:
        print(f"Species {species} is not available.")
        raise

    def extract_value(s, key):
        for item in s:
            if key in item:
                return item.split()[1].replace('"', "")
        return None

    url = (
        f"https://duffel.rail.bio/recount3/{species}/annotations/gene_sums/"
        f"{species}.gene_sums.{version}.gtf.gz"
    )
    gene_anno = pd.read_csv(
        url,
        delimiter="\t",
        skiprows=5,
        header=None,
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )

    # First split the attribute column into a list
    gene_anno["split_attribute"] = gene_anno["attribute"].str.split(";")

    # Apply the extract_value function for each column
    gene_anno["gene_id"] = gene_anno["split_attribute"].apply(lambda x: extract_value(x, "gene_id"))
    gene_anno["gene_type"] = gene_anno["split_attribute"].apply(
        lambda x: extract_value(x, "gene_type")
    )
    gene_anno["gene_name"] = gene_anno["split_attribute"].apply(
        lambda x: extract_value(x, "gene_name")
    )

    # Drop the unnecessary columns
    gene_anno = gene_anno.drop(["attribute", "split_attribute"], axis=1)
    return gene_anno


def get_cache_info(cache_type: str):
    """Get information about the cache

    Parameters
    ----------
    cache_type : str
        The entity of cached metadata, either 'samples', 'project' or 'genes'

    Returns
    -------
    str
        The usage of the current cache for the chosen entity.

    Raises
    ------
    ValueError
        If cache_type is not supported

    """
    if cache_type == "samples":
        return download_sample_metadata.cache_info()
    elif cache_type == "project":
        return download_project_metadata.cache_info()
    elif cache_type == "genes":
        return download_gene_metadata.cache_info()
    else:
        raise ValueError(f"Cache for {cache_type} is not available.")
