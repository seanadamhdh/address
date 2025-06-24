ADDRESS
================
Anita Alexandra Sanchez, Sean Adam
2025-06-24

- [IMPORTANT MESSAGE (UPDATE)](#important-message-update)
- [ADDRESS data exchange and online
  repository](#address-data-exchange-and-online-repository)
  - [GENERAL UPLOAD RULES /
    GUIDELINES](#general-upload-rules--guidelines)
  - [Spectral processing and model
    calibration](#spectral-processing-and-model-calibration)
  - [Using the models](#using-the-models)
- [Literature](#literature)

# IMPORTANT MESSAGE (UPDATE)

(New) For better access **moved the ADDRESS project folder** that
contains all data relevant for modelling **to Biogeochemie zfs1
server**. zfs1/…/hydropedo/projects/ADDRESS -\>
zfs1/…/fak3biogeochemie/03 Projects - Projekte/ADDRESS/ADDRESS:

In zfs1/…/hydropedo/projects/ADDRESS, there is now a shortcut linking to
the new location for easier access from the soil team side.

- /GIS -\> zfs1/…/fak3biogeochemie/projects/ADDRESS/

- /data -\> zfs1/…/fak3biogeochemie/projects/ADDRESS/

- /models -\> zfs1/…/fak3biogeochemie/projects/ADDRESS/

- /R_main

  - -\> gitlab.hrz.tu-chemnitz.de/sa56kotu-at-tu-freiberg-de/ADDRESS (no
    longer maintained, I will not have access after my university
    contract concludes)

  - UPDATE: migrated to -\>
    gitlab.hrz.tu-chemnitz.de/seanadamhdh/address as public repo  

# ADDRESS data exchange and online repository

This repository serves as a platform for code exchange and its
collaborative developement for the ADDRESS project.

## GENERAL UPLOAD RULES / GUIDELINES

Please do not directly upload to /main. Generally stable additions and
clean data can be added to /work_in_progress. Unstable scripts /messy
data should be commited to seperate (new) branch first to avoid
contamination of main repository. Folder structure:

- R_main  
  Contails all scripts. Clean / working scripts stored directly. This is
  the only folder that is pushed to github.

  - temp  
    Working output directory. Plots, model outputs etc. should be saved
    here initially. This folder is .gitignore and will therefore only
    exist locally.

  - old imported scripts  
    For saving / dumping uncleaned scripts that were not directly used
    for data processing but are still related.

— Moved to f./fak3/biogeochemie/projects/ADDRESS/

- models  
  Clean, chekced and evaluated model outputs.

- data  
  Raw data and compiled datasets

  - raw  
    Store raw data from the lab or the spectrolyzers here.

  - processed  
    Consolidated datasets and other data products.

## Spectral processing and model calibration

Chemometric models were calibrated from UV-Vis scans
(“/data/processed/allspecoriginal_Oct2023.xlsx”) and corresponding
anlytical data (“/data/Autosampler_A12_clean.csv”) to predict various
properties, i.e., DOC and trace element concentrations, from in situ
Spectrolyzer data.

### Spectral processing

UV-Vis spectra were recorded from 200 to 720 nm at a resolution of 0.5
nm. Raw spectra were cleaned by removing implausible scans (absorbance
\< 0 or absorbance \> 4.5). Spectra were processed using the prospectr
package Stevens and Ramirez-Lopez (2022). First, spectra were smoothed
using a Savitzky Golay filter with a polynomial degree of 3 and a window
size of 11 to remove noise. Then, a Standard Normal Variate (SNV)
correction was applied. Both raw, smoothed and smoothed + snv-corrected
spectra were used for model calibration.

### Model calibration

Chemometric models were calibrated using the Cubist package Kuhn and
Quinlan (2023) in the caret framework Kuhn (2008) for the following
variables: EC, pH, doc_mgL, dic_mgL, Sr, Suva254, Mn_mgL, Ni_mgL,
As_mgL, Cd_mgL, Pb_mgL, Fe_mgL, Al_mgL, Cu_mgL, Zn_mgL, and Durchfluss.
Models were calibrated both for untransformed and log1p-transformed
variables, and using the three differently processed spectral sets (raw,
smoothed, smoothed+snv). Cubist models were tuned for 1, 2, 5, 10, 20,
and 50 committees and 0-9 neighbours using a 10-fold cross-validation.
75 % of the dataset was used for calibration and 25 % was held back for
testing. The dataset was split using random, percentile-binned sampling.
Model accuracy was evaluated using the test set, using a modified
version of the evaluate_model() function found in the simplerspec
package Baumann (2020) to calculate RMSE, R2, RPD, and Lin’s concordance
correlation coefficient among other valuation statistics. An example for
the predictions for the test set can be seen in the figure below for
DOC.

![](images/log1p_spc_DOC.svg)

Full evaluation summary: “//zfs1.hrz.tu-freiberg.de/fak3biogeochemie/03
Projects - Projekte/ADDRESS/ADDRESS/models/Cubist_2024_02_21_eval.csv”

## Using the models

In the repository, ./R_main/load_specData_highRes.R contains a function
`predict_spectrolyzer()` for loading the collected spectrolyzer data and
applying the calibrated models. THe function returns a vector with
predictions. Each entry corresponds to a row of the input X. Example
use:

``` r
doc_mgL_pred=predict_spectrolyzer(X=spectrolyzer_data$spc_sg, 
                    model_dir="//zfs1.hrz.tu-freiberg.de/fak3biogeochemie/03 Projects - Projekte/ADDRESS/ADDRESS/models/Cubist_2024-02-21/",
                    variable="DOC_mgL", #select variable to predict
                    trans="log1p", # select transformation (predictions are already transformed back)
                    set="spc_sg", # must be the the same as set selcted as `X`
                    prefix="cubist-auto"  # default for model run
)

# i.e., add to dataset
spectrolyzer_data$DOC_predictions=doc_mgL_pred
```

# Literature

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-baumann2020" class="csl-entry">

Baumann, Philipp. 2020. “Simplerspec: Soil and Plant Spectroscopic Model
Building and Prediction.”
<https://github.com/philipp-baumann/simplerspec>.

</div>

<div id="ref-kuhn2008" class="csl-entry">

Kuhn, Max. 2008. “Building Predictive Models in r Using the Caret
Package.” *Journal of Statistical Software* 28 (5).
<https://doi.org/10.18637/jss.v028.i05>.

</div>

<div id="ref-kuhn2023" class="csl-entry">

Kuhn, Max, and Ross Quinlan. 2023. “Cubist: Rule- and Instance-Based
Regression Modeling.” <https://CRAN.R-project.org/package=Cubist>.

</div>

<div id="ref-stevens2022" class="csl-entry">

Stevens, Antoine, and Leonardo Ramirez-Lopez. 2022. “An Introduction to
the Prospectr Package.”
<https://cran.r-project.org/web/packages/prospectr/index.html>.

</div>

</div>
