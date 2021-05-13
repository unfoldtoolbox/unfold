# UNFOLD - TOOLBOX

A toolbox for deconvolution of overlapping EEG signals and (non)-linear modeling

* Linear deconvolution
* Model specification using R-style formulas (EEG~1+face+age)
* Programmed in a modular fashion
* Spline regression 
* Regularization (using glmnet)
* Temporal basis functions (Fourier & Splines)
* Estimate temporal response functions (TRFs) for time-continuous predictors
* Cross-validation

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
Check out the [toolbox tutorials](https://www.unfoldtoolbox.org/toolboxtutorials.html) for more information!
```
EEG = tutorial_simulate_data('2x2')
EEG = uf_designmat(EEG,'eventtypes',{'fixation'},'formula','y ~ 1+ cat(stimulusType)*cat(color)')
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.5 1])
EEG = uf_glmfit(EEG)
% (strictly speaking optional, but recommended)
ufresult = uf_condense(EEG)
ax = uf_plotParam(ufresult,'channel',1);
```

### Citation
Please cite as:

Ehinger BV, Dimigen O: "Unfold: An integrated toolbox for overlap correction, non-linear modeling, and regression-based EEG analysis",  peerJ 2019, https://doi.org/10.7717/peerj.7838
