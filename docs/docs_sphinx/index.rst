.. unfold documentation master file
.. highlight:: matlab

.. container:: landing

  .. image:: ../../media/unfold_800x377.png
    :align: center
    :width: 50 %

  A toolbox for *deconvolution* of overlapping EEG signals and *(non)-linear modeling*

  * Overlap between subsequent potentials using linear deconvolution
  * Model specification using R-style formulas, e.g. ``EEG~1+face+age``
  * Non-linear effects using regression splines, e.g. ``EEG~1+face+spl(age,10)``
  * Multiple events, e.g. *Stimulus*, *Response* and *Fixation*
  * Programmed in a modular fashion
  * Temporal basis functions (Fourier & Splines)
  * Optional regularization using **glmnet**


.. toctree::
   :maxdepth: 4
   :hidden:

   overview
   contenttutorials
   toolboxtutorials
   documentation
