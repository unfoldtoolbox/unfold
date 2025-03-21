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

> [!NOTE]
>  Have a look at the [Unfold.jl Julia ecosystem](https://github.com/unfoldtoolbox/UnfoldDocs/) which includes all the latest features, bugfixes and more. Unfold-Matlab is not actively developed, only bugfix releases are provided. PRs for new features are definitely welcome!
> 
### Getting help
  📢 **Try out our [discussion forum](https://github.com/unfoldtoolbox/unfold/discussions)** - we often get questions via email, a more transparent and open way is to use the [github discussions feature](https://github.com/unfoldtoolbox/unfold/discussions)


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

In addition, consider also citing the following reference, which illustrates the possibilites and options of unfold for a specific application example:

Dimigen O, Ehinger BV: "Regression-based analysis of combined EEG and eye-tracking data: Theory and applications. Journal of Vision, 21(1), 3-3",  https://jov.arvojournals.org/article.aspx?articleid=2772164


### Research notice
Please note that this repository is participating in a study into sustainability
 of open source projects. Data will be gathered about this repository for
 approximately the next 12 months, starting from June 2021.

Data collected will include number of contributors, number of PRs, time taken to
 close/merge these PRs, and issues closed.

For more information, please visit
[our informational page](https://sustainable-open-science-and-software.github.io/) or download our [participant information sheet](https://sustainable-open-science-and-software.github.io/assets/PIS_sustainable_software.pdf).
