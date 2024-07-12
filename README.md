This respository contains all of the scripts needed to perform the analysis and create the plots from [arXiv:2205.13546](https://arxiv.org/abs/2205.13546). 

To use this package, run [`setup.py`](setup.py) with the ``develop`` or
``install`` option.

The [`aidm/`](aidm/) directory contains the main calculation scripts. Adding
your own atom interferometer experiments can be done in the
[`aidm/experiments.py`](aidm/experiments.py) script.

The [`drivers/'](drivers/) directory contains the scripts needed to obtain the
limits and make the plots seen in the paper.

NOTE: The code does allow you to run decoherence rates for the atom
interferometer and BEC experiments. However, these rates wrongly include an N^2
enhancement that does not exist for these experiments. We recommend only
running the phase measurements for these types of experiments.


