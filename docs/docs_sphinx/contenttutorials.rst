Content Tutorials
==================
We recommend the articles of Smith on the regression-based framework of ERP waveforms (`Link Smith & Kutas 2015a <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5308234/>`,b) for an introduction to linear modeling in ERP analyses. The second part explicitly discusses deconvolution. This overview will just give you the gist on these topics.

Introduction to linear modeling
----------------------------------
Linear modeling tries to explain the data (ERPs / EEG signals) as a sum of some predictors / independent variables. In M/EEG analysis usually one GLM is fitted separately for each electrode and time-point. It can be formally stated as

.. math::
  y = x_1*\beta_1 + x_2*\beta_2 + ... + x_n*\beta_n + e = X\beta +e

because  :math:`\beta` is unknown one has to estimate it, usually using:

.. math::
  \beta = X^{-1}y


The choice of predictors (:math:`x_1,x_n...`) is up to the analyst of the signal. In order to understand how to use the toolbox it is helpful to understand the designmatrix :math:`X`. This designmatrix has multiple columns, each representing a predictor. Each row is one repetition, usually in EEG, one trial. If the predictor is a categorical one, some kind of encoding has to be used (reference/effects are currently supported). This encoding translates for example a two condition predictor 'face' vs 'house' in one predictor called the intercept, which in reference coding represents the reference-category (we use 'face') and the other predictor represents the difference of 'face' and 'house', because everything 'face' can explain, has already been soaked up by the first predictor.

Introduction to additive/spline linear modeling
------------------------------------------------
Sometimes relations between predictors and lets say P100 are not linear. For example saccadic amplitude has a logarithmic relationship. If this relationship is known, one can simply transform the predictor (log transform in this case) and then perform a linear fit. But in many cases either the relationship is not known, or it is not a simple function.

.. image:: ../../docs/tutorials/nonlinear1.png

In this figure we see that a linear function does not fit well our logarithmical relation.

In this case additive linear modeling allows to use flexible basis functions to model those non-linearities. A simple basis function is the boxcar function and it has a simple analogy: It is equivalent to splitting a continuous predictor in multiple categorical ones.

.. image:: ../../docs/tutorials/nonlinear_boxcar.png

In this example we split the continuous predictor into seven categorical predictors. In the modelfit (right plot) one can clearly see the step-function of this approach. A more sensible approach would enforce smooth borders and try to control the smoothness of the function. This can be achieved by using spline-basis functions instead of boxcars

.. image:: ../../docs/tutorials/nonlinear_spline.png

In order to get from the basisfunctions (left) to the function fit (right), each basis-function is multiplied by a fitted beta-coefficient value and then summed. These weighted basis-functions are in addition plotted in the right plot. It is important to note that the number of basis-functions is important to prevent over or underfitting. In the unfold toolbox one has to set the number of splines by hand. Nested crossvalidation to get a good estimate of the number of splines to use is certainly possible but computationally extremely expensive. In the field of additive modeling this issue is so far an unresolved problem.

Why deconvolution?
--------------------
It is standard practice in M/EEG research to separate events in time, so that no overlap occurs. But in some cases this is not possible. Especially in reading or free-viewing, that is experiments involving eye movements. Fixations occur on average every 250ms. Long lasting potentials like P300 or N400 occur around 400ms after stimulus onset. It is clear that in each fixation, the activity of the previous fixation is overlapping. Also even without long lasting potentials, baselines are influenced by the previous stimulus (e.g. Lütkenhof 2010).

If the overlap for two conditions is exactly identical, the differences between the conditions are unaffected. But often this is not the case. For example fixation durations are correlated with many cognitive variables of interest. And if fixation durations differ between two conditions, bias can arise as the overlap of the potentials is different. This can hide as well as generate differences between conditions.

How does deconvolution work?
----------------------------
Deconvolution tries to disentangle the signal. This is possible as the overlap is a bit different in each trial. Thus an algorithm can try to disentangle the process. The simplest way to do this is linear deconvolution, assuming that the signals of two events are linearly added. This has support in physiology: voltages are added linearly, that is, if two sources are active at the same time their voltage potentials can be summed. Linear deconvolution assumes that the first event does not influence the processing of the second event, an assumption that is most certainly wrong in the general case. But, importantly, one has to deal with this in a more generic ERP analysis as well. At least with deconvolution, some part of the overlap can be corrected for.

XX Explain timeshifting
