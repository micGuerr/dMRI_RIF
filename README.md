# Rotation Invariant Features for dMRI

## Purpose

In this repository we present a MATLAB implementation to obtain rotation invariant features from dMRI signal. The implementation is based on the [Zucchelli et al. (2019)](https://doi.org/10.1016/j.media.2019.101597) paper.

## Problem

dMRI is a powerful technique to infer microstructural properties of the living tissues in-vivo and non-invasively. The dMRI signal can be acquired with many different experimental settings which sensitize the signal to different properties of the tissues. For example, different gradient directions sensitize the signal to the microstructural environment along that specific direction.
Once such rich set of data has been acquired, it can then be used to infer the tissue microstructural properties, for example, via complicated biophysical models.

While in many applications this directional dependence is desirable (e.g. tractography), it may be not needed in others as it introduces signal changes which are independent of the quantity that we want to estimate (e.g. neurite density estimation).

Several techniques have been developed to remove the directional dependency of the diffusion signal and improve microstructural indexes estimation which take advantage of a rotationally invariant representation of the diffusion signal to factor out its directional dependency.

In this repository we present a MATLAB implementation of Zucchelli et al. framework for analytically generating a complete set of algebraically independent rotation invariant features (RIF) of the dMRI signal, given its Laplace-series expansion.

## Prerequisite

The framework has been developed and tested on MATLAB 2021a. The following packages and toolboxes are needed:

* The MATLAB [constrained spherical deconvolution](https://github.com/jdtournier/csd) toolbox.

If you want to compute the Gaunt matrices on your own, you'll need to have the:

* [Real/Complex Spherical Harmonic Transform, Gaunt Coefficients and Rotations](https://uk.mathworks.com/matlabcentral/fileexchange/43856-real-complex-spherical-harmonic-transform-gaunt-coefficients-and-rotations?s_tid=prof_contriblnk) from Aalto university.

If you want to compute the set of algebraically independent RIF on your own you'll need to have the:

* [MATLAB symbolic computation toolbox](https://uk.mathworks.com/help/symbolic/symbolic-computations-in-matlab.html).

## Download

You can download a zip version of the repository or use the following command from terminal:

```
git clone https://github.com/micGuerr/dMRI_RIF.git
```

## Installation

Once downloaded or cloned, you should navigate to the repository. Next, you should copy the `dmriRif_config.txt` file and save it with a `.m` extension.
Finallly, you should fill the scripts following the explenations contained in the file itself.

## Usage

### Input

### Output

### Example Command
