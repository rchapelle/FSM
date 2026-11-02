# Fit score method for sequential data synthesis

This repository provides the source code of the experiments for the paper: [A goodness of fit approach to select the optimal visit sequence for nonparametric data synthesis](https://arxiv.org/abs/2610.38559).

## Install the R package

Our method is available within the R package [opencesp](https://cran.r-project.org/web/packages/opencesp/index.html) available on CRAN. Easily install it with:  

```
install.packages("opencesp")
```
The actual experiments were performed with an early, development version of the package. To ensure that the results of the paper can be strictly replicated, we provide this early version of the opencesp package in the `opencesp.zip` archive. This version of the package can be built using the R development tools.

## Experiments

To run the experiments, first clone the repository with:
```
git clone https://github.com/rchapelle/FSM.git
```
Then, install the required R packages listed at the beginning of the `script.R` file, e.g.:
```
install.packages("randomForestSRC")
```
Finally, in the `script.R` file, set the variable `wd` to the directory in which the repository was cloned. All experiments are contained in the `script.R` file and can be conducted by running the script, for example within RStudio.

## Datasets

The [WBCD](https://doi.org/10.24432/C5HP4Z) dataset is licensed under the Creative Commons Attribution 4.0 International license.

## Reference

If you use this repository in your research, please cite our work using the following reference:

```
@misc{chapelle2026goodnessfitapproachselect,
      title={A goodness of fit approach to select the optimal visit sequence for nonparametric data synthesis}, 
      author={R{é}my Chapelle and Bruno Falissard and Mohammed Sedki},
      year={2026},
      eprint={2610.38559},
      archivePrefix={arXiv},
      primaryClass={stat.ML},
      url={https://arxiv.org/abs/2610.38559}, 
}
```
