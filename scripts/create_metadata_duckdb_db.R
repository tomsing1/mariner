library(dplyr)
library(here)
library(duckdb)
library(stringr)

db_file <- here::here("..", "data", "recount3.db")

# retrieve recount3 metadata from github
projects <- local({
  temp_file <- tempfile(fileext = ".Rdata")
  url = paste0("https://github.com/LieberInstitute/recount3-docs/raw/",
               "refs/heads/master/study-explorer/projects_meta.Rdata")
  download.file(url, temp_file)
  new_env <- new.env()
  load(temp_file, envir = new_env)
  get("projects_meta", envir = new_env)
}) |>
  dplyr::mutate(n_samples = as.integer(n_samples)) |>
  dplyr::select("project", "study_title", "study_abstract") |>
  dplyr::distinct()

count_files <- read.csv(
  paste0("https://github.com/LieberInstitute/recount3-docs/raw/refs/heads/",
         "master/study-explorer/",
         "recount3_raw_project_files_with_default_annotation.csv")
) |>
  dplyr::mutate(
    annotation = stringr::str_match(
      string = gene, 
      pattern = "gene_sums\\.\\w*\\.(\\w*)\\.gz$")[, 2]) |>
  dplyr::select("project", "organism","annotation", "file_source", 
                "n_samples", "gene")

metadata_files <- read.csv(
  paste0("https://github.com/LieberInstitute/recount3-docs/raw/refs/heads/",
         "master/study-explorer/recount3_metadata_files.csv")
) |> 
  dplyr::mutate(
    project = stringr::str_match(string = project_meta, 
                                 pattern = "metadata/\\w*/(\\w*)/")[, 2],
    file_source = dplyr::case_when(
      grepl("/gtex", project_meta, fixed = TRUE) ~ "gtex",
      grepl("/tcga", project_meta, fixed = TRUE) ~ "tcga",
      .default = "sra"
    ),
    organism = stringr::str_match(
      string = project_meta, 
      pattern = "recount3/(\\w*)/data_sources")[, 2]) |>
  dplyr::select("project", "file_source", "organism", "project_meta", 
                "recount_project")

files <- dplyr::inner_join(
  count_files, 
  metadata_files,
  by = join_by(project, organism, file_source)
)

# create SQLite database
con <- dbConnect(duckdb::duckdb(), db_file)

# add project table
dbExecute(
  con, 
  paste(
    "CREATE TABLE tbl_project", 
    "(", 
    "project TEXT PRIMARY KEY,",
    "study_title TEXT,",
    "study_abstract TEXT",
    ")")
)
dbAppendTable(con, name = "tbl_project", value = projects)


# add file table (some projects have files from both mouse and human)
dbExecute(
  con, 
  paste(
    "CREATE TABLE tbl_file", 
    "(", 
    "project TEXT,",
    "organism TEXT,",
    "annotation TEXT,",
    "file_source TEXT,",
    "n_samples INTEGER,",
    "recount_project TEXT,",
    "gene TEXT,",  
    "project_type TEXT,",
    "project_meta TEXT,",
    "FOREIGN KEY(project) REFERENCES tbl_project(project)",
    ")")
)
dbAppendTable(
  con, 
  name = "tbl_file", 
  value = files)
dbDisconnect(con)


# sandbox
if (FALSE ) {
  con <- dbConnect(RSQLite::SQLite(), db_file)
  dbListTables(con)

  # fully joined table
  dplyr::inner_join(tbl(con, "tbl_project"), tbl(con, "tbl_file"),
                    by = join_by("project")) |> 
    dplyr::collect() 

  # projects with more than one file (= human and mouse data available)  
  dplyr::inner_join(tbl(con, "tbl_project"), tbl(con, "tbl_file"),
                    by = join_by("project")) |>
    dplyr::group_by(project) |>
    dplyr::count() |>
    dplyr::filter(n > 1) |>
    dplyr::collect()
  
  dbDisconnect(con)
}