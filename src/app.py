import marimo

__generated_with = "0.11.13"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import utils
    return mo, pd, utils


@app.cell
def _(utils):
    samples = utils.download_sample_metadata("SRP115307")
    samples
    return (samples,)


@app.cell
def _(samples):
    samples.sample_title
    return


@app.cell
def _(samples):
    samples.sample_attributes
    return


@app.cell
def _(pd, samples):
    def split_attributes(attr_string):
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
                attributes[key.strip()] = value.strip()

        return attributes


    # Apply the function to create a new dataframe with split attributes
    attributes_df = pd.DataFrame(
        [split_attributes(row) for row in samples["sample_attributes"]],
        index=samples.index,
    )

    # Drop columns with constant values
    single_value_cols = [
        col for col in attributes_df.columns if attributes_df[col].nunique() == 1
    ]
    attributes_df.drop(columns=single_value_cols, inplace=True)
    samples_expanded = pd.concat(
        [samples[["experiment_acc", "sample_title"]], attributes_df], axis=1
    )
    samples_expanded
    return attributes_df, samples_expanded, single_value_cols, split_attributes


@app.cell
def _(utils):
    project = utils.download_project_metadata("SRP115307")
    project
    return (project,)


@app.cell
def _(mo):
    mo.md(r"""The project metadata contains the `organism` column, which we can used to retrieve the correct gene annotation (for mouse or human genes).""")
    return


@app.cell
def _(project):
    organism = project.organism.unique()
    assert len(organism) == 1
    organism = organism[0]
    organism
    return (organism,)


@app.cell
def _(organism, utils):
    genes = utils.download_gene_metadata(organism=organism)
    genes = genes.set_index("gene_id", drop=False)
    genes.head()
    return (genes,)


@app.cell
def _(utils):
    counts = utils.download_counts("SRP115307")
    counts.head()
    return (counts,)


@app.cell
def _():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    from anndata import AnnData
    import numpy as np
    return AnnData, DefaultInference, DeseqDataSet, DeseqStats, np


@app.cell
def _(counts, genes, samples_expanded):
    metadata = samples_expanded.loc[counts.columns]
    metadata["condition"] = metadata.gender
    keep_samples = ~metadata.condition.isna()
    metadata = metadata.loc[keep_samples]

    annotated_genes = list(set(counts.index).intersection(genes.index))
    genes2 = genes.loc[annotated_genes]

    m = counts.loc[annotated_genes].to_numpy().T
    m = m[keep_samples.tolist(), :]
    return annotated_genes, genes2, keep_samples, m, metadata


@app.cell
def _(DefaultInference, DeseqDataSet, genes, m, metadata):
    inference = DefaultInference(n_cpus=4)
    dds = DeseqDataSet(
        counts=m,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
    )
    dds.vars = genes
    dds.deseq2()
    dds
    return dds, inference


@app.cell
def _(DeseqStats, dds, inference):
    ds = DeseqStats(
        dds, contrast=["condition", "male", "female"], inference=inference
    )
    ds.summary()
    return (ds,)


@app.cell
def _(ds, pd):
    df = pd.concat([ds.dds.vars.reset_index(drop=True), 
                    ds.results_df.reset_index()], axis=1)
    df = df[~df.pvalue.isna()]
    df = df.sort_values("pvalue", ascending=True)
    df
    return (df,)


@app.cell
def _(ds):
    ds.dds.vars.shape
    return


if __name__ == "__main__":
    app.run()
