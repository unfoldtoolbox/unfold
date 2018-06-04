Unfold's toolbox workflow
=========================

We start with continuous data with events that mark times in the data that are of interest. Usually these are stimulus onsets, buttonpresses etc.

A minimal example is the following

.. code-block:: matlab

  EEG = uf_designmat('eventtypes',{'fixation'},'formula','1+cat(category)')
  EEG = uf_timeexpandDesignmat('timelimits',[-0.5 1])
  EEG = uf_glmfit(EEG)
  % (strictly speaking optional, but recommended)
  ufresult = uf_condense(EEG)



Definition of the design (:func:`uf_designmat`)
  We first need to define which event should be included in the deconvolution. Let's assume we have an event type 'stimulus' and one 'keypress'. We first define a formula for each event. These formulas are commonly used in `R` and also matlab and are called Wilkinson-Formulas. They allow to easily specify model-designs, e.g. 2x2 categorical linear model with interaction; continuous variables, spline-based (additive models, GAM) and mixtures between these models. The unfold-toolbox takes care of generating the appropriate designmatrix (denoted as `X`) for your design-specifications.


Timeexpansion of the designmatrix (:func:`uf_timeexpandDesignmat`)
  In order to perform the deconvolution, we need to expand the designmatrix over time. Because these designmatrices can become quite large, we offer the option to use a basis set (either cubic-splines or fourier-components). This will greatly reduce the number of parameters that need to be estimated but is an advanced feature - the benefits/costs are not well explored. The resulting timeexpanded/expanded designmatrix is denoted as `dcX`.


Fitting the model (:func:`uf_glmfit`)
  There are multiple ways to fit the expanded designmatrix `dcX`. The unfold-toolbox offers a memory-friendly iterative least squares algorithm (lsmr), default matlab pseudo-inverse and regularized glmnet (LASSO,ridge-regression and elastic-net).


Condensing Results (:func:`uf_condense`)
  This step extracts the betas, transforms them if necessary, and returns a cleaned-up output structure. When using time-basis functions (splines/fourier) we need to transform the estimated betas back to their native time-space. For example if we use 10 time-splines instead of 100 sample points, we receive 10 parameter-estimate (betas). In this case we want to transform the 10 betas back to the domain of 100 samples.


Plotting (e.g. :func:`uf_plotParam`)
  The toolbox offers multiple functions to explore also quite complex designs. These functions work on 1D (ERP-like), 2D (erpimage-like) or as timeseries of topoplots.

Group-Level Statistics
  Most commonly users estimate the parameters for all subjects and then test the parameters against an H0 of no effect. If no prior knowledege on the time/location of the effect is known, we recommend to use cluster-permutation based statistics. We recommend the EPT-TFCE toolbox. We focus on single subject beta-estimates and leave the statistics up to the user.

Minimal Data Specifications
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


Bad & Missing data
=========================================
There are some important differences in the workflow to a classical EEG experiment. This mostly regards to the handling of bad data (artefacts) and the removal/estimation of missing predictor values of some events (imputation)

Removing artifactual data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The toolbox requires continuous, uncleaned data. Removing data before the deconvolution would possibly also remove events that overlap partially in 'clean' events. Therefore the classical removal of continuous data, or removal of epochs does not work directly.

The function :func:`uf_continuousArtifactExclude` allows one to reject continuous portions of data that are were marked on the continuous signal. We expect a 'winrej'-matrix, the same format as the matrices used by eeglab (that is: columns sample start, sample end and each row one segment). The function then removes (puts to 0) the entries of EEG.unfold.Xdc corresponding to these time points. They are thereby effectively removed from being modeled.


Imputation of missing data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In linear models, missing data need to be imputed ('interpolated') or the event needs to be excluded. We currently offer four methods to deal with this in the function :func:`uf_imputeMissing` (to be called after the design specification):

drop
  this removes events that have a missing field in a single explanatory variable

marginal
  fills in a random sample from the values in other events. This preserves the overall distribution of missing events

mean
  fill in missing predictor-values with the mean of the remaining events.

median
  fill in missing predictor-values with the median of the remaining events.


Massive univariate modeling (rERP)
=================================================================
It is possible to use the toolbox for massive univariate linear modeling (also known as regression ERP, rERP) without any deconvolution applied. A minimal script looks like this:

.. code-block:: matlab

  EEG = uf_designmat('eventtypes',{'fixation'},'formula','1+cat(category)')
  EEG = uf_epoch(EEG,'timelimits',[-0.5 1])
  EEG = uf_glmfit_nodc(EEG)
  % (strictly speaking optional, but recommended)
  ufresult = uf_condense(EEG)


In addition, we offer several functions to compare deconvolved and non-deconvolved analyses. See tutorial :doc:`toolbox-tut04`.
