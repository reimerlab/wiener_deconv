## weiner_deconv
### Wiener deconvolution methods for estimating concentration dynamics from fluorescence traces

### Files:
#### Noncausal Wiener deconvolution (deconv.py)
These are functions for computing Wiener deconvolution and denoised reconstruction of a 1D (e.g. fluorescence) time-series as in (1). Here, concentrations are estimated using a time-domain representation of the
impulse response for the fluorophore (see text in 1). The noise spectrum is assumed to be flat.

This code also includes a method for estimating the regularization parameter W based on desired absolute error. 

#### Causal Wiener deconvolution (causal_deconv.py)
This code was not used in (1) but is included here for research purposes. Causal Wiener deconvolution was implemented following the derivation in (2):

#### Test notebook (test_deconv.py)
Notebook for testing the noncausal deconvolution method

#### Fluorescence data (fluorescence_detrended.csv)
Detrended fluorescence data for use in the test notebook

1) Neyhart, E et al. (2024; accepted) Cortical acetylcholine dynamics are predicted by cholinergic axon activity and behavior state. Cell Reports.
2) Lindell, Elisabeth. (2004) "On the Influence of Sensor Dynamics in Estimation of Reactor Noise Signals." 
