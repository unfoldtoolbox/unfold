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
git clone https://github.com/unfoldtoolbox/unfold
git submodule update --init --recursive --remote
```

### Running
```
run('init_unfold.m')
```


### Simple example
Check out the [toolbox tutorials](https://www.unfoldtoolbox.org/docs_sphinx/_build/html/toolboxtutorials.html) for more information!
```
EEG = tutorial_simulate_data('2x2')
EEG = uf_designmat('eventtypes',{'fixation'},'formula','y ~ 1+ cat(stimulusType)*cat(color)')
EEG = uf_timeexpandDesignmat('timelimits',[-0.5 1])
EEG = uf_glmfit(EEG)
% (strictly speaking optional, but recommended)
ufresult = uf_condense(EEG)
ax = uf_plotParam(ufresult,'channel',1);
```

### Citation
Please cite as:

Ehinger BV, Dimigen O: "Unfold: An integrated toolbox for overlap correction, non-linear modeling, and regression-based EEG analysis", https://doi.org/10.1101/360156
