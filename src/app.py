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
    from functools import lru_cache
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    from pydeseq2.preprocessing import deseq2_norm
    from sqlalchemy import create_engine
    return (
        AnnData,
        DefaultInference,
        DeseqDataSet,
        DeseqStats,
        alt,
        create_engine,
        deseq2_norm,
        lru_cache,
        mo,
        np,
        pd,
        utils,
    )


@app.cell
def _(create_engine):
    DATABASE_URL = "sqlite:///data/recount3.sqlite"
    engine = create_engine(DATABASE_URL)
    return DATABASE_URL, engine


@app.cell
def _(mo):
    mo.md(
        r"""
        ## The recount3 project

        [recount3](https://rna.recount.bio/) provides access to uniformly processed RNA-seq data from more than 750,000 human and mouse samples. The raw data was retrieved from public repositories, e.g.NCBI's short read archive (SRA), the Genotype-Tissue Expression (GTEx) project and The Cancer Genome Atlas (TCGA). Wilks et al described their full workflow in their [2021 Genome Research open access paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02533-6).

        The datasets can be explored through a [custom shiny web application](https://jhubiostatistics.shinyapps.io/recount3-study-explorer/), and the companion [recount Bioconductor package](https://www.bioconductor.org/packages/release/bioc/html/recount.html) facilitates import of the data into R sessions for downstream analysis.

        The processed data, including gene-level, exon-level and exon-junction counts, are available for download as CSV files, along with sample metadata available from the original repositories. 

        To explore the available datasets, I have retrieved the recount3 metadata from its [github repository](https://github.com/LieberInstitute/recount3-docs/tree/master/study-explorer) and collated them in a simple SQLite database. Let's get a list of all human or mouse datasets, along with their titles and abstracts.

        💡 Click on the column names to sort & filter the content, and select a row to see details about each study.
        """
    )
    return


@app.cell
def _(mo):
    species_radio = mo.ui.radio(
        options=["Human", "Mouse"], value="Human", inline=True
    )
    max_samples = mo.ui.range_slider(
        start=2, stop=1000, value=[10, 200], label="Number of samples"
    )
    mo.hstack([species_radio, max_samples])
    return max_samples, species_radio


@app.cell
def _(engine, max_samples, mo, species_radio, tbl_file, tbl_project):
    projects = mo.sql(
        f"""
        SELECT 
            prj.project AS project, 
            study_title AS title, 
            study_abstract AS abstract,
            fil.organism AS species,
            fil.n_samples AS n,
            fil.gene AS gene,
            fil.project_meta AS metadata, 
            fil.recount_project AS samples
        FROM tbl_project AS prj 
        LEFT JOIN tbl_file AS fil ON prj.project = fil.project 
        WHERE fil.organism  = '{species_radio.value.lower()}'
        AND fil.n_samples BETWEEN {max_samples.value[0]} AND {max_samples.value[1]}
        ORDER BY n;
        """,
        engine=engine,
        output=False,
    )
    project_table = mo.ui.table(
        projects[["project", "n", "title", "abstract"]],
        selection="single",
        initial_selection=[0],
    )
    project_table
    return project_table, projects


@app.cell
def _(mo, tbl_file, tbl_project):
    from dataclasses import dataclass
    from sqlalchemy.engine.base import Engine


    @dataclass
    class Project:
        """Class for keeping track of project metadata."""

        project: str
        title: str
        abstract: str
        species: str
        n: int
        genes: str
        samples: str
        meta: str


    def create_project(project: str, species: str, engine: Engine):
        project_tuple = mo.sql(
            f"""
            SELECT 
                prj.project AS project, 
                study_title AS title, 
                study_abstract AS abstract,
                fil.organism AS species,
                fil.n_samples AS n,
                fil.gene AS gene,
                fil.project_meta AS metadata, 
                fil.recount_project AS samples
            FROM tbl_project AS prj 
            LEFT JOIN tbl_file AS fil ON prj.project = fil.project 
            WHERE fil.project  = '{project}'
            AND fil.organism  = '{species}'
            LIMIT 1;
            """,
            engine=engine,
            output=False,
        ).row(0)
        return Project(*project_tuple)
    return Engine, Project, create_project, dataclass


@app.cell
def _(create_project, engine, mo, project_table, species_radio):
    if project_table.value["project"][0]:
        proj = create_project(
            project=project_table.value["project"][0],
            species=species_radio.value.lower(),
            engine=engine,
        )
        summary = mo.md(f"""
        ## {proj.project}: {proj.title}

        - Species: {proj.species}
        - {proj.n} samples

        {proj.abstract}
        """)
        run_button = mo.ui.run_button(label="Retrieve data")
    return proj, run_button, summary


@app.cell
def _(mo, run_button, summary):
    mo.vstack([summary, run_button])
    return


@app.cell
def _(mo):
    mo.md(r"""To retrieve the gene-level counts and the available sample & gene annotations for the selected study, please click the `Retrieve data` button above.""")
    return


@app.cell
def _(mo, proj, run_button, utils):
    mo.stop(not run_button.value)
    with mo.status.spinner(
        subtitle=f"Downloading sample information ..."
    ) as _spinner:
        samples = utils.download_sample_metadata(proj.samples)
    with mo.status.spinner(
        subtitle=f"Downloading gene annotations ..."
    ) as _spinner:
        genes = utils.download_gene_metadata(proj.genes)
    with mo.status.spinner(
        subtitle=f"Downloading gene-level counts ..."
    ) as _spinner:
        counts = utils.download_counts(proj.genes)
    return counts, genes, samples


@app.cell
def _(mo, proj, samples):
    sample_fields = [x for x in samples.columns if not x in ("experiment_acc")]
    mo.md(f"""
    The {proj.n} samples in this project are annotated with the following metadata fields:\n\n {"\n\n".join(["- " + item for item in sample_fields])}
    """)
    return (sample_fields,)


@app.cell
def _(sample_fields, samples):
    samples[sample_fields]
    return


@app.cell
def _(mo):
    mo.md(r"""For differential expression analysis, please choose the covariate of interest. Optionally, you may also select one or more covariates to adjust for, e.g. these variables will be included in the linear model along with the condition of interset.""")
    return


@app.cell
def _(mo, samples):
    group = mo.ui.dropdown(
        options=samples.columns, label="Choose condition of interest"
    )
    group
    return (group,)


@app.cell
def _(group, mo, samples):
    adjust_for = mo.ui.dropdown(
        options=[x for x in samples.columns if x != group.value],
        label="Choose covariates to adjust for",
    )
    adjust_for
    return (adjust_for,)


@app.cell
def _(adjust_for, group):
    if group.value:
        design = f"~{group.value}"
        if adjust_for.value:
            design = " + ".join([f"~{group.value}", adjust_for.value])
    return (design,)


@app.cell
def _(design, mo):
    mo.md(
        f"""Great! Please click the button below to fit a DESeq2 model to the dataset, 
        using your chosen design:  `{design}`"""
    )
    return


@app.cell
def _(AnnData, counts, genes, group, mo, samples):
    def create_annData(counts, samples, genes, column):
        # subset samples by excluding missing values in the `gender` column
        metadata = samples.loc[counts.columns]
        keep_samples = ~metadata[column].isna()
        metadata = metadata.loc[keep_samples]

        # subset genes to those that are annotated as 'protein_coding'
        annotated_genes = list(set(counts.index).intersection(genes.index))
        genes = genes.loc[annotated_genes]
        genes = genes.loc[genes.gene_type == "protein_coding"]

        # coerce the counts data.frame to a numpy array for faster transposition
        m = counts.loc[genes.index].to_numpy().T
        m = m[keep_samples.tolist(), :]
        adata = AnnData(X=m, obs=metadata, var=genes)
        return adata


    if group.value:
        with mo.status.spinner(subtitle=f"Created annData object ...") as _spinner:
            adata = create_annData(counts, samples, genes, group.value)
    return adata, create_annData


@app.cell
def _(mo):
    fit_button = mo.ui.run_button(label="Fit DESeq2 model")
    fit_button
    return (fit_button,)


@app.cell
def _(group):
    contrast = group.value
    return (contrast,)


@app.cell
def _(contrast, mo):
    mo.md(
        f"""Next, please select the two categorical levels of the `{contrast}` condition you would like to compare:"""
    )
    return


@app.cell
def _(adata, group, mo):
    level_1 = mo.ui.dropdown(
        options=adata.obs[group.value].unique(),
        label="Numerator",
    )
    return (level_1,)


@app.cell
def _(adata, group, level_1, mo):
    level_2 = mo.ui.dropdown(
        options=[x for x in adata.obs[group.value].unique() if x != level_1.value],
        label="Denominator",
    )
    level_selection = mo.ui.batch(
        mo.md(
            """
        - Numerator: {level_1}
        - Denominator: {level_2}
        """
        ),
        {"level_1": level_1, "level_2": level_2},
    ).form()
    level_selection
    return level_2, level_selection


@app.cell
def _(DefaultInference, DeseqDataSet, adata, design, fit_button, group, mo):
    mo.stop(not fit_button.value)


    def fit_model(adata, design, n_cpus=4, refit_cooks=True):
        inference = DefaultInference(n_cpus=n_cpus)
        dds = DeseqDataSet(
            adata=adata, design=design, refit_cooks=True, inference=inference
        )
        dds.deseq2()
        return dds


    if group.value:
        with mo.status.spinner(
            subtitle=f"Fitting DESeq2 model with design {design} ..."
        ) as _spinner:
            dds = fit_model(adata, design=design)
    return dds, fit_model


@app.cell
def _(DeseqStats, contrast, dds, level_selection, mo):
    def _(contrast, levels):
        with mo.status.spinner(
            subtitle=f"Extracting Wald test p-values ..."
        ) as _spinner:
            ds = DeseqStats(
                dds,
                contrast=[contrast] + levels,
                cooks_filter=True,
            )
            ds.summary()  # creates ds.results.df
            return ds


    if level_selection.value:
        chosen_levels = list(level_selection.value.values())
        ds = _(contrast, chosen_levels)
    return chosen_levels, ds


@app.cell
def _(chosen_levels, contrast, mo):
    mo.md(f"""
    The following table shows the differential expression results comparing 
    `{contrast}`:`{chosen_levels[-1]}` with `{contrast}`:`{chosen_levels[-2]}`. 
    """)
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
        initial_selection=[0],
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
def _(mo):
    mo.md(r"""To visualize the normalized or raw counts for one or more genes, please select one or more rows in the result table above, and press `Submit`.""")
    return


@app.cell
def _(dds, group, mo):
    covariate = mo.ui.dropdown(
        options=dds.obs_keys(),
        value=group.value,
        label="",
    )
    normalized = mo.ui.radio(
        options=["Normalized", "Raw"], value="Normalized", inline=True
    )
    markdown = mo.md(
        """
        - x-axis: {covariate}
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
        if batch.value:
            if batch.value["count_type"] == "Normalized":
                adata.X = deseq2_norm(adata.X)[0]
            ps = [
                scatter_plot(
                    adata=adata, gene=g, covariate=batch.value["covariate"]
                )
                for g in stat_table.value["gene_name"]
            ]
            return facet_wrap(ps, 3)


    _()
    return


if __name__ == "__main__":
    app.run()
