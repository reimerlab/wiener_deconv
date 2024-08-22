#### ```wiener_deconv```: Wiener deconvolution methods for estimating concentration dynamics from fluorescence traces

##### Noncausal Wiener deconvolution (```deconv.py```)
These are functions for computing Wiener deconvolution and denoised reconstruction of a 1D (e.g. fluorescence) time-series as in ref. 1. Here, concentrations are estimated using a time-domain representation of the
impulse response for the fluorophore. The noise spectrum is assumed to be flat.

This code also includes a method for estimating the regularization parameter W based on desired absolute error. 

##### Test notebook (```test_deconv.py```)
Notebook for testing the noncausal deconvolution method.

##### Fluorescence data (```fluorescence_detrended.csv```)
Detrended fluorescence data for use in the test notebook.

##### Causal Wiener deconvolution (```causal_deconv.py```)
Causal Wiener deconvolution was also implemented following the derivation in ref. 2. These methods were not used in ref. 1 but are included here for research purposes. 

#### References
1) Neyhart, E et al. (2024; accepted) Cortical acetylcholine dynamics are predicted by cholinergic axon activity and behavior state. Cell Reports.
2) Lindell, Elisabeth. (2004) "On the Influence of Sensor Dynamics in Estimation of Reactor Noise Signals." 
