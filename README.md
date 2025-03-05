# spmr
"spmr" is a R package for Causal Estimation and Inference in Nonlinear Mendelian Randomization Studies. When the causal effect of an exposure on an outcome is nonlinear, traditional MR estimation methods that rely on linear assumptions become unsuitable. To address this issue, we developed spmr.

If the functional form of the causality is known in advance, users can use two-stage prediction method or control function method to estimate the coefficients and standard errors of the causal function, thereby enabling causal inference. When the functional form of the causality is not well understood, semiparametric estimation can be employed to directly estimate the shape of the causal effect from the data, followed by the corresponding causal inference. Furthermore, if the instrumental variables exhibit horizontal pleiotropy, users can still obtain consistent estimates by adjusting the relevant parameters in "spmr".

For statistical details, please refer to the article "Causal Estimation and Inference in Nonlinear Mendelian Randomization Studies" by Xinpei Wang et al. 


## Installation
```r
devtools::install_github("xinpeiw/spmr")
```

## Functions
| Function name | Description   | 
|--------|---------|
| `format_sp` | Format data for the two-stage prediction method | 
| `format_cf` | Format data for the control function method | 
| `format_cf_semipara` | Format data for the semi-parametric estimation method (spMR) | 
| `sp_para` | Conduct the two-stage prediction estimation | 
| `cf_para`  |  Conduct the control function estimation  | 
| `cf_semipara` | Conduct semi-parametric estimation (spMR) | 

type `?function_name` to get more detailed parameter information for each function.

#### Feel free to reach out to Xinpei (xinpeiw at uchicago dot edu) if you have any problems :)
