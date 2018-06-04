Data Structures
====================================


Minimal Data Specifications
----------------------------
Data need to be in continuous form. We use the default format of EEG.

The minimal required fields (the usual EEGlab-structure) are:

* EEG.data (chan x time)

* EEG.times (vector of time points in milliseconds)

* EEG.srate (samplingrate in Hz)

* EEG.event (structure of events, same structure as in EEGlab)

  * event.type (e.g. 'stimulus' or 'keypress')
  * event.latency (in samples)
  * event.customField (`customField` = name of your variable, e.g. `color`. The content should be in this case a string 'red' or, for continuous variables, a number)




Fields of `EEG.unfold`
----------------------

splines
  all information on the splines are saved in here (see below). Each spline is added at splines{end+1}

formula
  contains all formulas specified in :func:`uf_designmat`

X
  The designmatrix. This can be used for 'classical' mass-univariate linear modeling (:func:`uf_epoch` and :func:`uf_glmfit_nodc`)

variabletypes
  Type of each variable/predictor, can be 'categorical', 'interaction', 'continuous' and 'spline'

variablenames
  Name of each variable/predictor without modifiers for level / spline modifier (e.g. factorA, sac_amplitude)

colnames:
  The name of each column of `X`. This field contains the modifier for level and spline (e.g. factorA_face or sac_amplitude_3.5)

cols2variablenames
  A list connecting the columns of `X` with the variables.

cols2eventtypes
  A list connecting the columns of `X` with possibly multiple events.

eventtypes
  The names of the events that are modeled. Only interesting if multiple different events were modeled.

Xdc
  Timeexpanded designmatrix [nsamples x (npredictors x ntimebasisfunctions)]. Output of `uf_timeexpandDesignmat`. If you need to modify this, have a look at `uf_designmat_addcol` to see which fields should be modified.

timebasis
  The basis-function of the timeexpand for the deconvolution. This matrix could be the identity matrix in case of "stick"/dirac-functions. Used only for splines/fourier time-basis functions

times: [1×20 double]
  A vector containing the time in seconds over what range the timeexpand occurred. This encodes the time of the resulting ERP

Xdc_terms2cols
  A list connecting the columns of `Xdc` with columns of `X`

beta_dc
  deconvolved betas. Output of `uf_glmfit`. This is the main outcome of this toolbox

beta_nodc
    non-deconvolved betas. This is a mass univariate fit where each timepoint and each electrode were fitted independently. Output of :func:`uf_glmfit_nodc`

channel
  for which channel the deconvolved betas have been calculated

effects_mean
    In case of effects coding contains the mean of the designmatrix columns


Fields of `ufresult`
----------------------
the ufresult structure is the output of :func:`uf_condense`. This function removes the time-splines if used and possibly evaluates splines at (automatically) specific quantiles.

unfold
  same as EEG.unfold

times
  same as ufresult.unfold.basistime, thus the epoch-time in s

chanlocs
  same as EEG.chanlocs

param
  a structure defining for each beta-value which event, what predictor, which variable-type and what the corresponding value is.

beta_nodc
  the betas without deconvolution [channel x time x predictors]

beta
  the betas with deconvolution [channel x time x predictors]

Fields of `unfold.splines`
------------------------------------
paramValues
  the parameter values of each event, e.g. for saccade amplitude: [1.3, 2.3, 6, 1.2 ...]

nSplines
  the number of splines used for modelling

knots
  the knot sequence. This is necessary to evaluate splines at a later point in time

splineFunction
  the function used to define the spline, could be a custom function.

removedSplineIdx
  The index of the spline which was removed during spline-generation. It is necessary to remove one spline in order to not have any collinearities. Depending on configuration either a middle or the first spline is removed.

X
  the entries of X times the spline (i.e. the subset of X)

name
  name of the spline

colnames
  column names that the spline will get in EEG.unfold.X
