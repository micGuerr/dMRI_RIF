# Rotation Invariant Features for dMRI

## Purpose

In this repository we present a MATLAB implementation to obtain rotation invariant features from dMRI signal. The implementation is based on the [Zucchelli et al. (2019) paper](https://doi.org/10.1016/j.media.2019.101597).

## Problem

dMRI is a powerful technique to infer microstructural properties of the living tissues in-vivo and non-invasively. The dMRI signal can be acquired with many different experimental settings which sensitize the signal to different properties of the tissues. For example, different gradient directions ...
Once such rich set of data has been acquired, it can then be used to infer the tissue microstructural properties, for example, via complicated biophysical models.

While in many applications this directional dependence is desirable (e.g. tractography), it can be a hindrance (not needed) to others as it introduces signal changes which are independent of the quantity that we want to estimate (e.g. neurite density estimation).

Several techniques have been developed to remove the directional dependency of the diffusion signal and improve microstructural indexes estimation which take advantage of a rotationally invariant representation of the diffusion signal to factor out its directional dependency.

In this repository we present a MATLAB implementation of Zucchelli et al. framework for analytically generating a complete set of algebraically independent rotation invariant features (RIF) of the dMRI signal, given its Laplace-series expansion.

## Prerequisite

## Download

## Installation

## Usage

### Input

### Output

### Example Command
