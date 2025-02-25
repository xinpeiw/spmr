# spmr
"smpr" is a package for Causal Estimation and Inference in Nonlinear Mendelian Randomization Studies. When the causal effect of an exposure on an outcome is nonlinear, traditional MR estimation methods that rely on linear assumptions become unsuitable. To address this issue, we developed "spmr".

If the functional form of the causality is known in advance, one can use two-stage prediction method or the control function method to estimate the coefficients and standard errors of the causal function, thereby enabling causal inference. When the functional form of the causality is not well understood, semiparametric estimation can be employed to directly infer the shape of the causal effect from the data, followed by the corresponding causal inference. Furthermore, if the instrumental variables exhibit horizontal pleiotropy (correlated or uncorrelated), users can still obtain consistent estimates by adjusting the relevant parameters in "spmr".

For statistical details, please refer to the article "Causal Estimation and Inference in Nonlinear Mendelian Randomization Studies" by Xinpei Wang. Feel free to reach out to Xinpei Wang if you have any problems:)



## Installation
```r
devtools::install_github("xinpeiw/spmr")
