Unfold's toolbox workflow
=========================

We start with continuous data with events that mark times in the data that are of interest. Usually these are stimulus onsets, buttonpresses etc.


Definition of the design (`dc_designmat`)
  We first need to define which event should be included in the deconvolution. Let's assume we have an event type 'stimulus' and one 'keypress'. We first define a formula for each event. These formulas are commonly used in `R` and also matlab and are called Wilkinson-Formulas. They allow to easily specify model-designs, e.g. 2x2 categorical linear model with interaction; continuous variables, spline-based (additive models, GAM) and mixtures between these models. The unfold-toolbox takes care of generating the appropriate designmatrix (denoted as `X`) for your design-specifications.


Timeexpansion of the designmatrix (`dc_timeexpandDesignmat`)
  In order to perform the deconvolution, we need to expand the designmatrix over time. Because these designmatrices can become quite large, we offer the option to use a basis set (either cubic-splines or fourier-components). This will greatly reduce the number of parameters that need to be estimated. The resulting timeexpanded/expanded designmatrix is denoted as `dcX`.


Fitting the model (`dc_glmfit`)
  There are multiple ways to fit the expanded designmatrix `dcX`. The unfold-toolbox offers a memory-friendly iterative least squares algorithm (lsmr), default matlab pseudo-inverse and regularized glmnet (LASSO,ridge-regression and elastic-net).


Extracting parameters (dc_beta2unfold)
  This step extracts the betas, transforms them if necessary, and returns a cleaned-up output structure. When using time-basis functions (splines/fourier) we need to transform the estimated betas back to their native time-space. For example if we use 10 time-splines instead of 100 sample points, we receive 10 parameter-estimate (betas). In this case we want to transform the 10 betas back to the domain of 100 samples.


Plotting (e.g. `dc_plotParam`)
  The toolbox offers multiple functions to explore also quite complex designs. These functions work on 1D (ERP-like), 2D (erpimage-like) or as timeseries of topoplots.

Group-Level Statistics
  Most commonly users estimate the parameters for all subjects and then test the parameters against an H0 of no effect. If no prior knowledege on the time/location of the effect is known, we recommend to use cluster-permutation based statistics. We recommend the EPT-TFCE toolbox. We focus on single subject beta-estimates and leave the statistics up to the user.

Specification for Input-Data
=============================
Data need to be in continuous form. We use the default format of EEG.

The minimal required fields (the usual EEGlab-structure) are:

* EEG.data (chan x time)

* EEG.times (vector of time points in milliseconds)

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



splines
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

cols2eventtype
  A list connecting the colums of `X` with the event. Multiple events can be modeled but they are not separable based on `X` alone.

eventtype
  The names of the events that are modeled. Only interesting if multiple different events were modeled.

Xdc
  Timeexpanded designmatrix [nsamples x (npredictors x ntimebasisfunctions)]. Output of `dc_timeexpandDesignmat`. If you need to modify this, have a look at `dc_designmat_addcol` to see which fields should be modified.

timebasis
  The basis-function of the timeexpand for the deconvolution. This matrix could be the identity matrix in case of "stick"/dirac-functions. Useful only for splines/fourier time-basis functions

times: [1×20 double]
  A vector containing the time in seconds over what range the timeexpand occured. This encodes the time of the resulting ERP

Xdc_terms2cols
  A list connecting the columns of `Xdc` with columns of `X`

beta_dc
  deconvolved betas. Output of `dc_glmfit`. This is the main outcome of this toolbox

beta_nodc
    non-deconvolved betas. This is a massive univariate fit where each timepoint and each electrode were fitted independently. Output of `dc_glmfit_nodc`

channel
  for which channel the deconvolved betas have been calculated

effects_mean
    In case of effects coding contains the mean of the designmatrix columns


Fields of `unfold`
----------------------
the unfold structure is the output of `dc_beta2unfold`. This function removes the time-splines if used and possibly evaluates splines at (automatically) specific quantiles.

deconv
  same as EEG.deconv

times
  same as EEG.times, thus the epoch-time in ms

chanlocs
  same as EEG.chanlocs

param
  a structure defining for each beta-value which event, what predictor, which variable-type and what the corresponding value is.

beta_nodc
  the betas without deconvolution [channel x time x predictors]

beta
  the betas with deconvolution [channel x time x predictors]

Fields of `deconv.splines`
------------------------------------
paramValues
  the parameter values of each event, e.g. for saccade amplitude: [1.3, 2.3, 6, 1.2 ...]
nSplines
  the number of splines used for modelling

knots
  the knot sequence. This is necessary to evaluate splines at a later point in time

removedSplineIdx
  The index of the spline which was removed during spline-generation. It is necessary to remove one spline in order to not have any collinearities. Depending on configuration either a middle or the first spline is removed.

X
  the entries of X times the spline (i.e. the subset of X)

name
  name of the spline

colnames
  column names that the spline will get in EEG.deconv.X
