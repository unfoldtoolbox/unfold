Unfold's toolbox workflow
=========================

We start with continuous data with events that mark times in the data that are of interest. Usually these are stimulus onsets, buttonpresses etc.


Definition of the design (`dc_designmat`)
  We first need to define which triggers (let's assume we have an event type 'stimulus' and one 'keypress') should be included in the deconvolution. Because each event could belong to a different condition, which again could be different for stimuli/keypresses, we can define a simple formula for each event. These formulas are commonly used in `R` and also matlab and are called Wilkinson-Formulas. They allow to easily specify very simple designs, e.g. 2x2 with only categorical, but also continuous variables, spline-based (additive models, GAM) and mixture between these models. The unfold-toolbox takes care of generating the appropriate designmatrix (denoted as `X`) for your design-specifications.


Timeshift of the designmatrix (`dc_timeexpandDesignmat`)
  In order to perform the deconvolution, we need to expand the designmatrix over time. Because these designmatrices can become quite large, we offer the option to use a basis set (either cubic-splines or fourier-components). This will greatly reduce the number of parameters that need to be estimated. The resulting timeshifted/expanded designmatrix is denoted as `dcX`.


Fitting the model (`dc_glmfit`)
  There are multiple ways to fit the expanded designmatrix `dcX`. The unfold-toolbox offers a memory-friendly iterative least squares algorithm, default matlab pseudo-inverse and regularized glmnet (LASSO,ridge-regression and elastic-net).


Transforming parameters (dc_beta2unfold)
  (optional) In some cases, e.g. using time-basis functions (splines/fourier) or when using spline-nonlinear predictors, it is helpful to transform the estimated parameters back to their native space. For example if we use 10 time-splines instead of 100 sample points, we receive 10 parameter-estimate (betas). In this case we want to transform the 10 betas back to the domain of 100 samples. This increases the interpretation. This step is also possible for continuous predictors, but because a slope is readily interpretable, not as necessary.


Plotting (e.g. `dc_plotParam`)
  The toolbox offers multiple functions to explore also quite complex designs. These functions work on 1D (ERP-like), 2D (erpimage-like) or as timeseries of topoplots.

Group-Level Statistics
  Most commonly users estimate the parameters for all subjects and then test the parameters against an H0 of no effect. If no prior knowledege on the time/location of the effect is known, we recommend to use cluster-permutation based statistics. We recommend the EPT-TFCE toolbox which we included in our toolbox for convenience. The toolbox currently does not have functions to do single-subject statistics.

Specification for Input-Data
=============================
Data need to be in continuous form. We use the default format of EEG.


The minimal required fields (the usual EEGlab-structure) are:

* EEG.data (chan x time)

* EEG.times (vector of time points)

* EEG.srate (samplingrate in Hz)

* EEG.event (structure of events, same structure as in EEGlab)

  * event.type (e.g. 'stimulus' or 'keypress')
  * event.latency (in samples)
  * event.customField (`customField` = name of your variable, e.g. `color`. The content should be in this case a string 'red' or, for continuous variables, a number)

Imputation of missing data
==================================
In linear models, missing data need to be imputed ('interpolated') or the event needs to be excluded. We currently offer four methods to deal with this in the function `dc_imputeMissing` (to be called after the design specification):

drop
  this removes events that have a missing field in a single explanatory variable

marginal
  fills in a random sample from the values in other events. This preserves the overall distribution of missing events

mean
  fill in missing predictor-values with the mean of the remaining events.

median
  fill in missing predictor-values with the median of the remaining events.


Removing artefactual data
============================
Because we start with continuous data, usual removal of data stretches does not work the same. The function `dc_continuousArtifactExclude` allows one to reject continuous portions of data that are either marked automatically or manually. We expect a 'winrej'-matrix, the same format as the matrices used by eeglab (that is: columns sample start, sample end and each row one segment)

Comparison between deconvolved and non-deconvolved (classical)
================================================================
we offer several functions to compare deconvolved and non-deconvolved analyses. See tutorial :doc:`toolbox-tut04`

Explanation of all variable fields
==================================
Fields of `EEG.deconv`
----------------------

effects_mean
  In case of effects coding contains the mean of the designmatrix columns. We do not want to keep the original event structure around, thus we save this information here

predictorSplines
  all information on the splines are saved in here (see below). Splines are a bit more involved to describe, additional information needs to be saved

formula
  contains all formulas specified in `dc_designmat`

X
  The designmatrix. This can be used for 'classical' massive-univariate linear modeling

variableType
  Type of each variable/predictor, can be 'categorical','continuous' and 'spline'

variableNames
  Name of each variable/predictor

colnames: {1×10 cell}
  The name of each column of `X`. As soon as categorical variables are used, we will have more columns than variables.

cols2variableNames
  A list connecting the columns of `X` with the variables. A 0 means that the intercept is meant.

col2eventtype
  A list connecting the colums of `X` with the event. Multiple events can be modeled but they are not separable based on `X` alone.

eventtype
  The names of the events that are modeled. Only interesting if multiple different events were modeled.

dcX
  Timeshifted designmatrix [nsamples x (npredictors x ntimebasisfunctions)]. Output of `dc_timeexpandDesignmat`. If you need to modify this, have a look at `dc_designmat_addcol` to see which fields should be modified.

dcBasis
  The basis-function of the timeshift for the deconvolution. This matrix could be the identity matrix in case of "stick"/dirac-functions. Useful only for splines/fourier time-basis functions

dcBasistime: [1×20 double]
  A vector containing the time in seconds over what range the timeshift occured. This encodes the time of the resulting ERP

dcX_termidx
  A list connecting the columns of `dcX` with columns of `X`

dcBeta
  deconvolved betas. Output of `dc_glmfit`. This is the main outcome of this toolbox

channel
  for which channel the deconvolved betas have been calculated

XBeta
  non-deconvolved betas. This is a massive univariate fit where each timepoint and each electrode were fitted independently. Output of `dc_glmfit_nodc`


Fields of `unfold`
----------------------
the unfold structure is the output of `dc_beta2unfold`. This function removes the time-splines if used and possibly evaluates splines at (automatically) specific quantiles.

deconv
  same as EEG.deconv

times
  same as EEG.times, thus the epoch-time in ms

chanlocs
  same as EEG.chanlocs

epoch
  a structure defining for each beta-value to what event and predictor it belong. In case 'convertsplines' was activated in `dc_beta2unfold` the spline-betas were evaluated at specified values, then these values are also saved here.

beta_nodc
  the betas without deconvolution [channel x time x predictors]
beta
  the betas with deconvolution [channel x time x predictors]

Fields of `deconv.predictorSplines`
------------------------------------
paramValues
  [5988×1 double]

nSplines
  9

min
  : -7.9975

max
  : 7.9782

spline2val
  : [1×1000 double]

spline2val_idx
  : [1×5988 double]

basis
  : [1000×9 double]

X
  : [5988×9 double]

name
  : 'splineA'

colnames
  : {1×9 cell}
