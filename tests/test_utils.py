from io import StringIO

import pandas as pd
import pytest

from src.utils import (
    Project,
    clean_column_names,
    download_counts,
    download_sample_metadata,
    extract_annotation_source,
    get_cache_info,
    sanitize,
    split_attributes,
)

# ---------- Basic unit tests ----------

@pytest.mark.parametrize("input_str,expected", [
    ("abc123", "abc123"),
    ("abc-123", "abc_123"),
    ("abc!@#123", "abc_123"),
    ("  abc   ", "abc"),
    ("---", ""),
])
def test_sanitize(input_str, expected):
    assert sanitize(input_str) == expected


def test_clean_column_names():
    df = pd.DataFrame(columns=["col 1", "col-2", "col@3"])
    cleaned = clean_column_names(df.copy())
    assert cleaned.columns.tolist() == ["col_1", "col_2", "col_3"]


@pytest.mark.parametrize("attr_string,expected", [
    ("age;;30|sex;;male", {"age": "30", "sex": "male"}),
    ("foo;;bar|baz;;qux", {"foo": "bar", "baz": "qux"}),
    ("no_separator", {}),
    (None, {}),
])
def test_split_attributes(attr_string, expected):
    result = split_attributes(attr_string)
    assert result == expected


@pytest.mark.parametrize("url,expected", [
    ("https://example.com/whatever.gene_sums.THING.M023.gz", "M023"),
    ("https://example.com/path.gene_sums.data.G026.gz", "G026"),
    ("no.matching.pattern", None)
])
def test_extract_annotation_source(url, expected):
    assert extract_annotation_source(url) == expected


def test_project_dataclass():
    project = Project(id="PRJ123", count_url="url1", sample_url="url2", project_url="url3")
    assert project.id == "PRJ123"
    assert project.count_url == "url1"
    assert project.sample_url == "url2"
    assert project.project_url == "url3"


# ---------- Mocking download functions ----------

@pytest.fixture
def fake_counts_file(monkeypatch):
    csv_content = "gene_id\tSRR1\tSRR2\nENSG1\t5\t10\nENSG2\t3\t7\n"
    real_read_csv = pd.read_csv
    def mock_read_csv(*args, **kwargs):
        return real_read_csv(StringIO(csv_content), delimiter="\t", index_col="gene_id")
    monkeypatch.setattr("src.utils.pd.read_csv", mock_read_csv)
    return mock_read_csv

def test_download_counts(fake_counts_file):
    df = download_counts("fake_url")
    assert df.shape == (2, 2)
    assert list(df.columns) == ["SRR1", "SRR2"]
    assert list(df.index) == ["ENSG1", "ENSG2"]


@pytest.fixture
def fake_metadata_file(monkeypatch):
    csv_content = "external_id\texperiment_acc\tsample_title\tsample_attributes"\
  "\nS1\tEXP1\ttitle1\tage;;30|sex;;male\nS2\tEXP2\ttitle2\tage;;31|sex;;female\n"
    real_read_csv = pd.read_csv
    def mock_read_csv(*args, **kwargs):
        return real_read_csv(StringIO(csv_content), delimiter="\t", index_col="external_id")
    monkeypatch.setattr("src.utils.pd.read_csv", mock_read_csv)
    return mock_read_csv

def test_download_sample_metadata(fake_metadata_file):
    df = download_sample_metadata("fake_url")
    assert "experiment_acc" in df.columns
    assert "sample_title" in df.columns
    assert "age" in df.columns
    assert "sex" in df.columns
    assert "age" in df.columns and df.shape[1] <= 4  # Drops single-value columns


def test_get_cache_info_valid():
    info = get_cache_info("samples")
    assert "hits" in str(info) or hasattr(info, "hits")  # Depending on Python version


def test_get_cache_info_invalid():
    with pytest.raises(ValueError):
        get_cache_info("not_a_valid_cache")


