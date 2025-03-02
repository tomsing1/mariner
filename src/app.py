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
    utils.download_sample_metadata("SRP115307")
    return


@app.cell
def _(utils):
    utils.download_project_metadata("SRP115307")
    return


@app.cell
def _(utils):
    utils.download_gene_metadata(species="mouse")
    return


@app.cell
def _(utils):
    utils.get_cache_info("samples")
    return


if __name__ == "__main__":
    app.run()
