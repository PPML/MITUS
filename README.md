# MITUS
Modeling Interventions for Tuberculosis in the United States 

### Project Goals
The goal of the MITUS package is to provide a transmission dynamic model to investigate the impacts of various TB control strategies in the United States. The model has been calibrated to national and state level data for a select number of states. Built-in control strategies are in three broad categories: 
- pre-defined scenarios
- programmatic changes to the LTBI and active TB diagnosis and care cascades
- targeted testing and treatment of customizable risk groups

### Model Description 

In the MITUS model, a core TB dimension captures TB transmission, natural history, and treatment. Additional dimensions represent 1) TB progression risk, 2) mortality risk, 3) socio-economic disadvantage, 4) LTBI treatment history, 5) nativity (U.S.-born or non-U.S.–born), and 6) age-based differences in disease mechanisms and risk factor prevalence. 
Below we include a schematic of the mathematical model of tuberculosis epidemiology in the United States. 

![MITUS state diagram](https://github.com/PPML/MITUS/blob/development/inst/MITUS_schematic.jpg)

### Key Data Sources 

**TB data**  

Reported TB incidence is collected from the publicly available data from the National Tuberculosis Surveillance System, which can be accessed through the [Online Tuberculosis Information System](https://wonder.cdc.gov/tb.html).
LTBI prevalence is based on an analysis of the 2011-2012 National Health and Nutrition Examination Survey by [Woodruff et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0140881) 

**Mortality data** 

Historical all cause mortality data are sourced from the [Human Mortality Database](https://mortality.org/Country/Country?cntr=USA)
Future projections of all cause mortality are sourced from the Social Security Administration 
Deaths with TB data are collected from National Center for Health Statistics [Multiple Cause of Death Data](https://wonder.cdc.gov/mcd.html)

**Population data** 

Population data is collected from the [American Community Survey](https://usa.ipums.org/usa/)

## Installation 
Currently, MITUS is available for download with an authentication key. In order to gain access to the package, please the author at [email](nswartwood@hsph.harvard.edu) with a detailed description of your interest in the package. 

```
# install.packages("devtools")
devtools::install_github("PPML/MITUS", auth_token='TOKEN')
```

