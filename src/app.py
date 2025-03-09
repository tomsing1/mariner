import marimo

__generated_with = "0.11.17"
app = marimo.App(width="medium")


@app.cell
def _():
    import altair as alt
    import marimo as mo
    import pandas as pd
    import numpy as np
    import utils

    from anndata import AnnData
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    from pydeseq2.preprocessing import deseq2_norm
    return (
        AnnData,
        DefaultInference,
        DeseqDataSet,
        DeseqStats,
        alt,
        deseq2_norm,
        mo,
        np,
        pd,
        utils,
    )


@app.cell
def _(utils):
    samples = utils.download_sample_metadata("SRP115307")
    return (samples,)


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
    return attributes_df, samples_expanded, single_value_cols, split_attributes


@app.cell
def _(utils):
    project = utils.download_project_metadata("SRP115307")
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
    return (organism,)


@app.cell
def _(organism, utils):
    genes = utils.download_gene_metadata(organism=organism)
    genes = genes.set_index("gene_id", drop=False)
    return (genes,)


@app.cell
def _(utils):
    counts = utils.download_counts("SRP115307")
    return (counts,)


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
    return (dds,)


@app.cell
def _(DeseqStats, dds):
    ds = DeseqStats(
        dds, contrast=["condition", "male", "female"], cooks_filter=True
    )
    ds.summary()  # creates ds.results.df
    return (ds,)


@app.cell
def _(mo):
    mo.md("""### Independent filtering""")
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
    stat_table = mo.ui.table(
        df.reset_index(drop=True),
        show_column_summaries=True,
        selection="multi",
        initial_selection=[1],
        label="Differential expression results",
    )
    stat_table
    return df, stat_table


@app.cell(hide_code=True)
def _(alt, pd):
    def facet_wrap(subplts, plots_per_row):
        rows = [
            subplts[i : i + plots_per_row]
            for i in range(0, len(subplts), plots_per_row)
        ]
        compound_chart = alt.hconcat()
        for r in rows:
            rowplot = alt.vconcat()  # start a new row
            for item in r:
                rowplot |= item  # add suplot to current row as a new column
            compound_chart &= rowplot  # add the entire row of plots as a new row
        return compound_chart


    def scatter_plot(adata, gene, covariate, jitter_scale=(-1, 2)):
        if gene in adata.var_names:
            gene_idx = adata.var_names.get_loc(gene)
        else:
            raise ValueError(f"Gene {gene} not found in var_names.")
        if not covariate in adata.obs_keys():
            raise ValueError(f"Covariate {covariate} not found in obs_keys.")

        _df = pd.DataFrame(
            {
                "normalized_counts": adata.X[:, gene_idx].flatten(),
                covariate: adata.obs[covariate].values,
            }
        )
        _df

        chart = (
            alt.Chart(_df)
            .mark_point(filled=True, size=100)
            .encode(
                x=alt.X(
                    f"{covariate}:N",
                    axis=alt.Axis(
                        title=covariate.title(),
                        labelFontSize=12,
                        titleFontSize=14,
                    ),
                ),
                y=alt.Y(
                    "normalized_counts:Q",
                    axis=alt.Axis(
                        title="Normalized counts",
                        labelFontSize=12,
                        titleFontSize=14,
                    ),
                ),
                xOffset=alt.XOffset(
                    "jitter:Q", scale=alt.Scale(domain=jitter_scale)
                ),
                color=alt.Color(f"{covariate}:N", legend=None),
            )
            .transform_calculate(jitter="random()")
            .resolve_scale()
            .properties(width=200, height=200)
            .properties(
                title=gene,
            )
            # .configure_title(fontSize=24)
        )
        return chart
    return facet_wrap, scatter_plot


@app.cell
def _(dds, mo):
    covariate = mo.ui.dropdown(
        options=dds.obs_keys(),
        value="condition" if "condition" in dds.obs_keys() else dds.obs_keys()[0],
        label="",
    )
    normalized = mo.ui.radio(
        options=["Normalized", "Raw"], value="Normalized", inline=True
    )
    markdown = mo.md(
        """
        - x-axis?: {covariate}
        - y-axis: {count_type}
        """
    )
    batch = mo.ui.batch(
        markdown, {"covariate": covariate, "count_type": normalized}
    ).form()
    batch
    return batch, covariate, markdown, normalized


@app.cell
def _(batch, dds, deseq2_norm, facet_wrap, scatter_plot, stat_table):
    def _():
        adata = dds.copy()
        adata.var = adata.var.set_index("gene_name", inplace=False, drop=False)
        if batch.value["count_type"] == "Normalized":
            adata.X = deseq2_norm(adata.X)[0]
        ps = [
            scatter_plot(adata=adata, gene=g, covariate=batch.value["covariate"])
            for g in stat_table.value["gene_name"]
        ]
        return facet_wrap(ps, 3)


    _()
    return


if __name__ == "__main__":
    app.run()
