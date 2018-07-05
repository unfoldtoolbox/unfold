|
|
.. unfold documentation master file
.. highlight:: matlab

.. container:: landing

  .. container:: center-div

    .. image:: ../../media/unfold_800x377.png
      :align: center
      :width: 30 %

Unfold - EEG Deconvolution Toolbox
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A toolbox for *deconvolution* of overlapping EEG signals and *(non)-linear modeling*

First draft of reference paper
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`Download our reference paper Ehinger & Dimigen 2018 from bioRxiv <https://www.biorxiv.org/content/early/2018/07/04/360156>`_

This is a draft version, we are very happy to receive comments and criticism (info@unfoldtoolbox.org)

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



.. toctree::
   :maxdepth: 4
   :hidden:

   overview
   contenttutorials
   toolboxtutorials
   documentation
   datastructures
   
