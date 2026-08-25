# The complex relationship between anti-immigrant sentiment and exposure in the Netherlands

This repository contains the code used to calculate spatial exposure to immigrants from CBS grid data and to reproduce the main data preparation, visualisation, and regression workflow for the accompanying paper.

The paper describes the research design, variables, models, and robustness analyses. This repository is intended to make the computational workflow inspectable and reusable, and to support replication attempts.

## Contents

* `exposure_fun.py` — functions for calculating distance-weighted exposure on CBS grid cells, aggregating cell-level measures to regions, merging election results, and preparing analytical variables.
* `Analysis_illustration.ipynb` — a Jupyter notebook illustrating the data-loading, exposure-calculation, aggregation, mapping, and regression workflow.
* `requirements.txt` — Python package versions used for this code.

## Setup

Create and activate a Python environment, then install the required packages:

```bash
python -m pip install -r requirements.txt
```

The code was developed with Python 3.11.

## Data

The analysis uses publicly available data from:

* **Statistics Netherlands (CBS):** 2023 500 m × 500 m grid data and 2023 administrative boundaries (`gebiedsindelingen`).
  * [https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data](https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data)
* **Kiesraad:** results of the 2023 Dutch House of Representatives election.
  * [https://www.verkiezingsuitslagen.nl/verkiezingen/detail/TK20231122](https://www.verkiezingsuitslagen.nl/verkiezingen/detail/TK20231122)

The paper also uses LISS Panel data to substantiate the link between voting behaviour and anti-immigrant sentiment. These data require a free LISS account:

* [https://www.dataarchive.lissdata.nl/study-units/view/22](https://www.dataarchive.lissdata.nl/study-units/view/22)

The main spatial analysis itself uses the CBS and Kiesraad data. Some election inputs in the notebook are pre-cleaned spatial or tabular files:

* `votes_stemB_2023.gpkg`
* `GM2021vote.gpkg`
* `wijk_votes_21.csv`

Place these files in the repository directory, or adjust the paths in the notebook. The two 2021 files are needed only for analyses of the change in right-wing vote share between 2021 and 2023.

The notebook expects the project folder and data folders to share a common parent directory:

```text
parent-folder/
├── Notebook-folder/
│   ├── Analysis_illustration.ipynb
│   ├── exposure_fun.py
│   ├── requirements.txt
│   ├── votes_stemB_2023.gpkg
│   ├── GM2021vote.gpkg
│   └── wijk_votes_21.csv
├── Grid_data/
│   ├── 2025-cbs_vk500_2023_v2/
│   │   └── cbs_vk500_2023_v2.gpkg
│   └── cbsgebiedsindelingen2016-2025/
│       └── cbsgebiedsindelingen2023.gpkg
├── Election_data/
│   └── Verkiezingsuitslag Tweede Kamer 2023 (CSV formaat)/
└── Gemeenten alfabetisch 2023.csv
```

From within `Notebook-folder/`, paths such as `../Grid_data/...` and `../Gemeenten alfabetisch 2023.csv` therefore refer to files and folders in the parent directory. You may instead organise the data differently and update the path variables near the beginning of the notebook.

## Running the workflow

Open `Analysis_illustration.ipynb` in Jupyter and run the cells in order.

1. Set the paths to the CBS grid data, CBS boundaries, election data, and municipality–province lookup file.
2. Select the unit of analysis: `wijk` (neighbourhood) or `gemeente` (municipality).
3. Set the exposure parameter `exposure_k`.
4. Run the exposure calculation, spatial aggregation, election-data merge, and desired regression or plotting sections.

For the paper’s main analyses, use:

```python
exposure_k = 10000
```

The notebook defaults to `1000` because this is quicker for illustration. The exposure calculation for the full country, particularly with a larger value of `k`, may take some time.

## Exposure measure

For each grid cell, the code expands through successive Manhattan-distance neighbourhood layers until the specified population threshold `k` is reached. It calculates a distance-weighted exposure measure using the proportion of residents born outside the Netherlands with a non-European heritage. The resulting cell-level measure is aggregated to neighbourhoods or municipalities, weighted by the Dutch population in each cell.

Further methodological detail is provided in the accompanying paper.

## Reuse and questions

We welcome replication attempts, questions, and suggestions for improving this repository. If you encounter an issue with data preparation, file paths, or interpretation of the workflow, please contact the corresponding author listed in the accompanying paper.
