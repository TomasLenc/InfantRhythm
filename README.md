Analysis code for the manuscript "_Infants show enhanced neural response to musical meter frequencies beyond low-level features_".

---

Dependencies:  
* [letswave6](https://github.com/NOCIONS/letswave6)
* [fieldtrip](https://www.fieldtriptoolbox.org/)

Data: 
* [source EEG data](https://gin.g-node.org/TomasLenc/XPInfant_source)
* [stimulus audio files](https://gin.g-node.org/TomasLenc/XPInfant_stimuli)

To run the whole pipeline, clone the current repository, rename the folder to `code` and place it in the project root directory. 

Download the EEG data from [GIN](https://gin.g-node.org/TomasLenc/XPInfant_source) (either using [datalad](https://www.datalad.org/) or simply by downloading the dataset as .zip). Extract the zip if needed, rename the folder to `source` and place into the root directory. 

Repeat the same procedure to download the [stimulus files](https://gin.g-node.org/TomasLenc/XPInfant_stimuli) and rename the folder to `stimuli`. 

Your project root directory should look like this: 
```
.
├── code
│   ├── lib
│   │   ├── AB
│   │   ├── AuditoryToolbox
│   │   ...
│   ├── chan2rm
│   ├── cochlear_model.m
│   ├── extractFeatures.m
│   ├── getParams.m
│   ├── main.m
│	...
│   ├── utils.R
│   └── XPInfant.Rproj
├── source
│   ├── high sync P001.lw6
│   ├── high sync P001.mat
│   ├── high sync P002.lw6
│   ...
├── stimuli
│   ├── H_syncopated_200ms_calib75dB_max_short-term_ff_freq1237Hz_onset10_offset50.wav 
│   ...
```

Make sure to update the paths in `getParams.m`. 

Then, run the whole analysis by executing the `main.m` script in matlab. 

---

Next step is to run statistical analyses in R. 

Open `XPInfant.Rproj` in RStudio and set your environment using the [renv](https://rstudio.github.io/renv/articles/collaborating.html) package. 

Update the path to your dataset root in `main.R` and run this script to generate reports and figures in the `derivatives/reports/` folder. 

