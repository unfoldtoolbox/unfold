|
|
.. unfold documentation master file
.. highlight:: matlab

.. container:: landing

  .. container:: center-div

    .. image:: ../../media/unfold_800x377.png
      :align: center
<<<<<<< HEAD
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
=======
      :width: 30 %



A toolbox for *deconvolution* of overlapping EEG signals and *(non)-linear modeling*


What can you do with *unfold*?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Adjust for **overlap** between subsequent potentials using linear deconvolution
* Massive-Univarite Modeling (rERP) using R-style formulas, e.g. ``EEG~1+face+age``
* Non-linear effects using **regression splines** (GAM), e.g. ``EEG~1+face+spl(age,10)``
* Model **multiple events**, e.g. *Stimulus*, *Response* and *Fixation*
* Use temporal basis functions (Fourier & Splines)
* (Optional) **regularization** using glmnet


Requirements
^^^^^^^^^^^^^^^^^
* MATLAB 2015a+
* Continuous data in EEGLAB 12+ format
* Unfold toolbox `Download it on GitHub <https://github.com/unfoldtoolbox/unfold/>`_

Getting Started
^^^^^^^^^^^^^^^^^
To get started, best is to start with the 2x2 ANOVA-Design tutorial :doc:`toolbox-tut01`


.. raw:: html
    <iframe style="border: 0; height: 200px; width: 600px;"       src="https://www.unfoldtoolbox.org/piwik/index.php?module=CoreAdminHome&action=optOut&language=en&backgroundColor=&fontColor=&fontSize=&fontFamily"></iframe>


>>>>>>> 8f4935e0e6877a23b6d0431490abfbf84a3f4d57


.. toctree::
   :maxdepth: 4
   :hidden:

   overview
   contenttutorials
   toolboxtutorials
   documentation
   datastructures
   impressum
