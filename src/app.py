import marimo

__generated_with = "0.11.13"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import utils

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    from anndata import AnnData
    import numpy as np
    return (
        AnnData,
        DefaultInference,
        DeseqDataSet,
        DeseqStats,
        mo,
        np,
        pd,
        utils,
    )


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
    return (genes,)


@app.cell
def _(utils):
    counts = utils.download_counts("SRP115307")
    counts.head()
    return (counts,)


@app.cell
def _():
    return


@app.cell
def _(
    AnnData,
    DefaultInference,
    DeseqDataSet,
    counts,
    genes,
    samples_expanded,
):
    def _(counts, samples, genes):
        # subset samples by excluding missing values in the `gender` column
        metadata = samples.loc[counts.columns]
        metadata["condition"] = metadata.gender
        keep_samples = ~metadata.condition.isna()
        metadata = metadata.loc[keep_samples]

        # subset genes to those that are annotated as 'protein_coding'
        annotated_genes = list(set(counts.index).intersection(genes.index))
        genes = genes.loc[annotated_genes]
        genes = genes.loc[genes.gene_type == "protein_coding"]

        # coerce the counts data.frame to a numpy array for faster transposition
        m = counts.loc[genes.index].to_numpy().T
        m = m[keep_samples.tolist(), :]

        inference = DefaultInference(n_cpus=4)
        dds = DeseqDataSet(
            adata=AnnData(X=m, obs=metadata, var=genes),
            design="~condition",
            refit_cooks=True,
            inference=inference,
        )
        dds.deseq2()
        return dds


    dds = _(counts=counts, samples=samples_expanded, genes=genes)
    dds
    return (dds,)


@app.cell
def _(DeseqStats, dds):
    ds = DeseqStats(
        dds, contrast=["condition", "male", "female"], cooks_filter=True
    )
    ds.summary()
    return (ds,)


@app.cell
def _(mo):
    mo.md(
        """
        ### Independent filtering


        """
    )
    return


@app.cell
def _(ds, mo, pd):
    df = pd.concat(
        [
            ds.dds.var[["gene_name", "seqname"]],
            ds.results_df,
        ],
        axis=1,
    )
    df = df[~df.pvalue.isna()]
    df = df.sort_values("pvalue", ascending=True)
    df = df.round(
        {
            "baseMean": 1,
            "log2FoldChange": 1,
            "lfcSE": 2,
            "stat": 1,
            "pvalue": 10,
            "padj": 10,
        }
    )
    mo.ui.table(df, show_column_summaries=True)
    return (df,)


if __name__ == "__main__":
    app.run()
