import pandas as pd
import re
from dataclasses import dataclass
from functools import lru_cache

GENE_ANNOTATIONS = {
  "M023": "mouse",
  "G026": "human"
}

@dataclass
class Project:
    """Class for keeping track of project metadata."""
    id: str
    count_url: str
    sample_url: str
    project_url: str


@lru_cache(maxsize=5)
def download_counts(url:str) -> pd.DataFrame:
    """
    Download and raw counts for a specific SRA project.

    Parameters
    ----------
    url : str
        The URL pointing to the recount3 gene counts CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing gene-level counts with
        one column for each run (SRR) in the dataset
    """

    return pd.read_csv(url, delimiter="\t", skiprows=2, index_col="gene_id")


@lru_cache(maxsize=10)
def download_metadata(url: str) -> pd.DataFrame:
    """
    Download and parse sample or project metadata for a specific SRA project.

    Parameters
    ----------
    url : str
        The URL pointing to the recount3 sample or project metadata CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing sample or project metadata
    """

    return pd.read_csv(url, delimiter="\t", index_col="external_id")


def sanitize(x: str) -> str:
  return re.sub('[^0-9a-zA-Z]+', '_', x.strip()).strip("_")


def clean_column_names(df):
    df.columns = [sanitize(x) for x in df.columns]
    return df

  
def split_attributes(attr_string: str):
    """Split a sample attributes string into a dictionary of key-value pairs.

    Args:
        attr_string (str): String containing attributes in format "key;;value|key;;value"

    Returns:
        dict: Dictionary of attribute key-value pairs
    """
    if pd.isna(attr_string):
        return {}

    pairs = attr_string.split("|")
    attributes = {}

    for pair in pairs:
        if ";;" in pair:
            key, value = pair.split(";;")
            attributes[sanitize(key)] = sanitize(value)

    return attributes


@lru_cache(maxsize=10)
def download_sample_metadata(url: str) -> pd.DataFrame:
    samples = download_metadata(url)
    attributes_df = pd.DataFrame(
        [split_attributes(row) for row in samples["sample_attributes"]],
        index=samples.index,
    )
    single_value_cols = [
        col for col in attributes_df.columns if attributes_df[col].nunique() == 1
    ]
    attributes_df.drop(columns=single_value_cols, inplace=True)
    df = pd.concat(
        [samples[["experiment_acc", "sample_title"]], attributes_df], axis=1
    )
    return clean_column_names(df)


def extract_annotation_source(url):
    """
    Extract the gene annotation identifier from the gene summary URL

    Parameters
    ----------
    url : str
        The URL pointing to the recount3 sample or project metadata CSV file.

    Returns
    -------
    str
        The gene annotation identifer (e.g. M023 or G026)
    """
    pattern = r"\.gene_sums\.[^\.]+\.([^.]+)\.gz"
    match = re.search(pattern, url)

    if match:
        return match.group(1)

  
@lru_cache(maxsize=10)
def download_gene_metadata(url) -> pd.DataFrame:
    """
    Download and parse gene annotations for a gene count file

    Parameters
    ----------
    url : str
        The URL pointing to the recount3 gene count file

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
    """
    annotation = extract_annotation_source(url)
    try:
        species = GENE_ANNOTATIONS[annotation]
    except KeyError:
        print(f"Annotation source {annotation} is not available.")
        raise

    url = (
        f"https://duffel.rail.bio/recount3/{species}/annotations/gene_sums/"
        f"{species}.gene_sums.{annotation}.gtf.gz"
    )
    def extract_value(s, key):
        for item in s:
            if key in item:
                return item.split()[1].replace('"', "")
        return None

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
        ]
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
    gene_anno.set_index('gene_id', drop=False, inplace=True)

    # retain only the first row for each gene symbol
    gene_anno.drop_duplicates(["gene_name"], inplace=True)
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
