# MITUS
Modeling Interventions of Tuberculosis in the United States 

### Project Goals
The goal of the MITUS package is to provide a deterministic model to investigate the impacts of various TB control strategies in the United States. The model has been calibrated to national and state level data for a select number of states. Built-in control strategies are in three broad categories: 
- pre-defined scenarios
- programmatic changes to the LTBI and active TB diagnosis and care cascades
- targeted testing and treatment of customizable risk groups

### Model Schematic 
Below we include a schematic of the deterministic mathematical model of tuberculosis epidemiology in the United States. 

![MITUS state diagram](https://raw.githubusercontent.com/PPML/MITUS/development/inst/MITUS_schematic.jpg?token=ABJHE33UJISBULG36F7GBG266ORCK)

## Installation 
Currently, MITUS is available for download with an authentication key. In order to gain access to the package, please [email](nswartwood@hsph.harvard.edu) with a detailed description of your interest in the package. 
```
# install.packages("devtools")
devtools::install_github("PPML/MITUS", auth_token='TOKEN')
```

## Code Examples 
For detailed code examples please see our vignettes. 
