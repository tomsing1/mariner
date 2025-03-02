import marimo

__generated_with = "0.11.13"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import utils
    return mo, pd, utils


@app.cell(disabled=True)
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
        [split_attributes(row) for row in samples["sample_attributes"]]
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
    utils.download_gene_metadata(organism=organism)
    return


@app.cell
def _(utils):
    utils.get_cache_info("samples")
    return


if __name__ == "__main__":
    app.run()
