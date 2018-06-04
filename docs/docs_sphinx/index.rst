.. unfold documentation master file
.. highlight:: matlab

.. container:: landing

  .. container:: center-div

    .. image:: ../../media/unfold_800x377.png
      :align: center
      :width: 25 %

  A Matlab toolbox for advanced regression-based EEG analysis *deconvolution* of overlapping EEG signals and *(non)-linear modeling*

What is the unfold toolbox?

Feature overview
  * Disentangle overlapping potentials from subsequent events using linear deconvolution
  * Control continuous linear or non-linear covariates
  * Easy model specification using R-style formulas, e.g. ``EEG~1+face+age``
  * Non-linear effects using regression splines, e.g. ``EEG~1+face+spl(age,10)``
  * Multiple events, e.g. *Stimulus*, *Response*, and *Fixation*
  * Programmed in a modular fashion
  * Use temporal basis functions (Fourier & Splines)
  * Optional regularization using **glmnet**

.. image:: ../tutorials/unfold_plotting_tools.png
  :align: center
  :width: 55 %


.. toctree::
   :maxdepth: 4
   :hidden:

   overview
   contenttutorials
   toolboxtutorials
   documentation
   impressum
