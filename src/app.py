import marimo

__generated_with = "0.11.31"
app = marimo.App(width="medium")


@app.cell
def _():
    from dataclasses import dataclass
    from functools import lru_cache

    import altair as alt
    import duckdb
    import marimo as mo
    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    from pydeseq2.preprocessing import deseq2_norm

    import utils
    return (
        AnnData,
        DefaultInference,
        DeseqDataSet,
        DeseqStats,
        alt,
        dataclass,
        deseq2_norm,
        duckdb,
        lru_cache,
        mo,
        np,
        pd,
        utils,
    )


@app.cell
def _():
    DATABASE_URL = (
        "https://github.com/tomsing1/mariner/raw/refs/heads/main/data/recount3.db"
    )
    return (DATABASE_URL,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## The recount3 project

        [recount3](https://rna.recount.bio/) provides access to uniformly processed RNA-seq data from more than 750,000 human and mouse samples. The raw data was retrieved from public repositories, e.g.NCBI's short read archive (SRA), the Genotype-Tissue Expression (GTEx) project and The Cancer Genome Atlas (TCGA). Wilks et al described their full workflow in their [2021 Genome Research open access paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02533-6).

        The datasets can be explored through a [custom shiny web application](https://jhubiostatistics.shinyapps.io/recount3-study-explorer/), and the companion [recount Bioconductor package](https://www.bioconductor.org/packages/release/bioc/html/recount.html) facilitates import of the data into R sessions for downstream analysis. The processed data, including gene-level, exon-level and exon-junction counts, are available for download as CSV files, along with sample metadata available from the original repositories. 

        To explore the available datasets, I have retrieved the recount3 metadata from its [github repository](https://github.com/LieberInstitute/recount3-docs/tree/master/study-explorer) and collated them in a 
        [duckdb database](https://github.com/tomsing1/mariner/raw/refs/heads/main/data/recount3.db). 
        Let's get a list of all human or mouse datasets, along with their titles and abstracts.

        💡 Click on the column names to sort & filter the content, and select a row to see details about each study.
        """  # noqa: E501
    )
    return


@app.cell
def _(duckdb, mo):
    con = duckdb.connect(database=":memory:")
    with mo.status.spinner(
        subtitle="Connecting to remote recount3 project database ..."
    ) as _spinner:
        con.sql("INSTALL httpfs;")
        con.sql(
            "ATTACH 'https://github.com/tomsing1/mariner/raw/refs/heads/main/data/recount3.db' "
            "as db;"
        )
    return (con,)


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
def _(con, db, max_samples, mo, species_radio, tbl_file, tbl_project):
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
            FROM db.tbl_project AS prj 
            LEFT JOIN db.tbl_file AS fil ON prj.project = fil.project 
            WHERE fil.organism  = '{species_radio.value.lower()}'
            AND fil.n_samples BETWEEN {max_samples.value[0]} AND {max_samples.value[1]}
            ORDER BY n;
            """,
        engine=con,
        output=False,
    )

    project_table = mo.ui.table(
        projects[["project", "n", "title", "abstract"]],
        selection="single",
        initial_selection=[0],
    )
    project_table  # noqa B018
    return project_table, projects


@app.cell
def _(dataclass, db, duckdb, mo, tbl_file, tbl_project):
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


    def create_project(
        project: str, species: str, engine: duckdb.duckdb.DuckDBPyConnection
    ):
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
            FROM db.tbl_project AS prj 
            LEFT JOIN db.tbl_file AS fil ON prj.project = fil.project 
            WHERE fil.project  = '{project}'
            AND fil.organism  = '{species}'
            LIMIT 1;
            """,
            engine=engine,
            output=False,
        ).row(0)
        return Project(*project_tuple)
    return Project, create_project


@app.cell
def _(con, create_project, mo, project_table, species_radio):
    if project_table.value["project"][0]:
        proj = create_project(
            project=project_table.value["project"][0],
            species=species_radio.value.lower(),
            engine=con,
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
    mo.md(
        r"""
        To retrieve the gene-level counts and the available sample & gene annotations for 
                  the selected study, please click the `Retrieve data` button above.
        """
    )
    return


@app.cell
def _(mo, proj, run_button, utils):
    mo.stop(not run_button.value)
    with mo.status.spinner(
        subtitle="Downloading sample information ..."
    ) as _spinner:
        samples = utils.download_sample_metadata(proj.samples)
    with mo.status.spinner(
        subtitle="Downloading gene annotations ..."
    ) as _spinner:
        genes = utils.download_gene_metadata(proj.genes)
    with mo.status.spinner(
        subtitle="Downloading gene-level counts ..."
    ) as _spinner:
        counts = utils.download_counts(proj.genes)
    return counts, genes, samples


@app.cell
def _(mo, proj, samples):
    sample_fields = [x for x in samples.columns if x not in ("experiment_acc")]
    sample_fields_string = "\n\n".join(["- " + item for item in sample_fields])
    mo.md(
        f"""The {proj.n} samples in this project are annotated with the following metadata fields:
        \n{sample_fields_string}"""
    )
    return sample_fields, sample_fields_string


@app.cell
def _(sample_fields, samples):
    samples[sample_fields]
    return


@app.cell
def _(mo, run_button, samples):
    mo.stop(not run_button.value)
    mo.md(
        "For differential expression analysis, please choose the covariate of interest. "
        "Optionally, you may also select one or more covariates to adjust for, e.g. "
        "these variables will be included in the linear model along with the condition of interest."
    )
    group = mo.ui.dropdown(
        options=samples.columns, label="Choose condition of interest"
    )
    group  # noqa B018
    return (group,)


@app.cell
def _(group, mo, samples):
    mo.stop(not group.value)
    adjust_for = mo.ui.dropdown(
        options=[x for x in samples.columns if x != group.value],
        label="Choose covariates to adjust for",
    )
    adjust_for  # noqa B018
    return (adjust_for,)


@app.cell
def _(AnnData, adjust_for, counts, genes, group, mo, samples):
    mo.stop(not group.value)
    contrast = group.value
    design = f"~{group.value}"
    if adjust_for.value:
        design = " + ".join([f"~{group.value}", adjust_for.value])


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


    mo.md(
        f"""Great! Please click the button below to fit a DESeq2 model to the dataset, 
        using your chosen design:  `{design}`"""
    )
    with mo.status.spinner(subtitle="Created annData object ...") as _spinner:
        adata = create_annData(counts, samples, genes, group.value)
    return adata, contrast, create_annData, design


@app.cell
def _(adata, group, mo):
    mo.stop(not group.value)
    level_1 = mo.ui.dropdown(
        options=adata.obs[group.value].unique(),
        label="Numerator",
    )
    return (level_1,)


@app.cell
def _(adata, group, level_1, mo):
    mo.stop(not group.value)
    level_2 = mo.ui.dropdown(
        options=[x for x in adata.obs[group.value].unique() if x != level_1.value],
        label="Denominator",
    )
    return (level_2,)


@app.cell
def _(level_1, level_2, mo):
    def validate_form(form):
        if not form["level_1"]:
            return "Please select the numerator of interest."
        if not form["level_2"]:
            return "Please select the denominator of interest."
        return None


    level_selection = mo.ui.batch(
        mo.md(
            """
        - {level_1}:
        - {level_2}:
        """
        ),
        {"level_1": level_1, "level_2": level_2},
    ).form(submit_button_label="Fit DESeq2 model", validate=validate_form)
    level_selection  # noqa B018
    return level_selection, validate_form


@app.cell
def _(DefaultInference, DeseqDataSet, adata, design, level_selection, mo, np):
    mo.stop(not level_selection.value)


    def fit_model(adata, design, n_cpus=4, refit_cooks=True):
        inference = DefaultInference(n_cpus=n_cpus)
        dds = DeseqDataSet(
            adata=adata, design=design, refit_cooks=True, inference=inference
        )
        try:
            dds.deseq2()
        except np.linalg.LinAlgError as err:
            if "Singular matrix" in str(err):
                return (
                    False,
                    "Sorry, that comparison is not possible because "
                    "the design matrix is singular. Please choose another set of covariates.",
                )
        return True, dds


    with mo.status.spinner(
        subtitle=f"Fitting DESeq2 model with design {design} ..."
    ) as _spinner:
        fit_success, dds = fit_model(adata, design=design)
    return dds, fit_model, fit_success


@app.cell
def _(dds, fit_success, mo):
    mo.md(
        f"{mo.icon('lucide:thumbs-up')} Successfully fit the DESeq2 model!"
    ) if fit_success else mo.md(f"{mo.icon('lucide:triangle-alert')}" + dds)
    return


@app.cell
def _(DeseqStats, contrast, dds, fit_success, level_selection, mo):
    mo.stop(
        not level_selection.value,
        mo.md("**Select the numerator and denomator to compare.**"),
    )

    mo.stop(
        not fit_success,
    )


    def _(contrast, levels):
        with mo.status.spinner(
            subtitle="Extracting Wald test p-values ..."
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
def _(chosen_levels, contrast, ds, fit_success, mo, pd):
    mo.stop(
        not fit_success,
    )
    mo.md(f"""
    The following table shows the differential expression results comparing 
    `{contrast}`:`{chosen_levels[-1]}` with `{contrast}`:`{chosen_levels[-2]}`. 
    """)
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
    stat_table  # noqa B018
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
        if covariate not in adata.obs_keys():
            raise ValueError(f"Covariate {covariate} not found in obs_keys.")

        _df = pd.DataFrame(
            {
                "normalized_counts": adata.X[:, gene_idx].flatten(),
                covariate: adata.obs[covariate].values,
            }
        )

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
def _(dds, fit_success, group, mo):
    mo.stop(
        not fit_success,
    )
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
        To visualize the normalized or raw counts for one or more genes, please select one or more
        rows in the result table above, and press `Submit`.

        Use the options below to customize the plots:

        - x-axis: {covariate}
        - y-axis: {count_type}
        """
    )
    batch = mo.ui.batch(
        markdown, {"covariate": covariate, "count_type": normalized}
    ).form()
    batch  # noqa B018
    return batch, covariate, markdown, normalized


@app.cell
def _(
    batch,
    dds,
    deseq2_norm,
    facet_wrap,
    fit_success,
    mo,
    scatter_plot,
    stat_table,
):
    mo.stop(
        not fit_success,
    )


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
