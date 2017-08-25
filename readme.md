# UNFOLD - TOOLBOX

A toolbox for deconvolution of overlapping EEG signals and (non)-linear modeling

* Linear deconvolution
* Model specification using R-style formulas (EEG~1+face+age)
* Programmed in a modular fashion
* Regularization using glmnet
* Spline-Regression
* Temporal basis functions (Fourier & Splines)

### Installation
``` 
git clone https://github.com/behinger/unfold
git submodule update --init --recursive --remote
```

### Running
```
run('init_unfold.m')
```
