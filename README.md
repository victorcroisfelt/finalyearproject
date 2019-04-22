# Spatial Correlation and Low Complexity Signal Processing Techniques in Massive MIMO Systems

This repository includes the PDF of the following final year project (FYP):

Victor Croisfelt Rodrigues, "Spatial Correlation and Low Complexity Signal Processing Techniques in Massive MIMO Systems", Final Year Project, Universidade Estadual de Londrina, Londrina, Brazil, December, 2018.

**Code:** This repository also has a research-oriented code package, located in the "matlab" folder, that allow readers to replicate the results of the technical report appended in Appendix A of the above FYP, entlited as:

Victor Croisfelt Rodrigues and Taufik Abrão, "Massive MIMO System in TDD Mode: Channel Estimation and Spectral Efficiency", Technical Report, Universidade Estadual de Londrina, Londrina, Brazil, December, 2018.

I hope this content helps in your reaseach and contributes to building the precepts behind *open science*. Remarkably, in order to boost the idea of open science and further drive the evolution of science, we also motivate you to share your published results to the public.

If you have any questions and if you have encountered any inconsistency, please do not hesitate to contact me via victorcroisfelt@gmail.com.

## FYP Abstract
One of the most promising technologies to match the current society requirements in the fifth generation (5G) of wireless communications is the massive multiple-input multiple-output (M-MIMO) system. At a first moment, this work focused to present the canonical concepts behind this technology in order to uncover its functionality and open issues related with its envisioned implementation. As a fundamental complication of this system, the pilot contamination was studied here; where a successive pilot decontamination method was duly assessed. Some drawbacks of this approach were pointed out, whereby the use of multiple pilot training phases was seen unattractive to handle with pilot contamination. Motivated by the lack of studies associated with spatial correlation over multi-antenna channels in an M-MIMO scenario, the exponential antenna array correlation model combined with large-scale fading variations over the array was analyzed for both channel hardening and favorable propagation effects; where, the normalized mean squared error (NMSE) was also deployed as an evaluation metric. Since the current interest resides on the use of different antenna array arrangements, the performance assessment of spatially correlated M-MIMO channels was carried out deploying the uniform linear array (ULA) and the uniform planar array (UPA). Generally speaking, the favorable propagation effect was found to improve with spatial correlation, wherein UPA guarantees better results. In view of the fact that favorable propagation is properly related to channel estimation, the pilot contamination was seen reduced when considering the use of the minimum mean-square error (MMSE) estimator. Instead of these positive results, channel hardening was observed to be poorly sustained in high-spatial correlation conditions. Given these points, a trade-off between the assurance of favorable propagation and channel hardening has been discussed in order to support a more realistic network design based on the studied model. Another current research interest lies in the reduction of the complexity associated with the linear signal processing techniques applicable to M-MIMO systems. With this in mind, the Kaczmarz algorithm (KA) was implemented to solve the combining/precoding problems following some recent directions found in the literature. So as to improve the rate of convergence of this iterative algorithm, a modification in the KA has been proposed, showing good results with respect to the original one.

## FYP Content
The final year project presented here is a composition of several results of which some of them have been published in the form of articles, we can cite the following:

- Victor Croisfelt Rodrigues and Taufik Abrão, "[An Evaluation of Successive Pilot Decontamination in Massive MIMO](http://www.uel.br/revistas/uel/index.php/semexatas/article/view/34450/24955)", Semina: Ciências Exatas e Tecnológicas, Londrina, Brazil, v. 39, n. 2, ago./dez. 2018. "[Code available here](https://github.com/victorcroisfelt/eval-suc-pilot-decon-mmimo)"   
- Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão, "[Exponential Spatial Correlation with Large-Scale Fading Variations in Massive MIMO Channel Estimation](https://doi.org/10.1002/ett.3563)", Transactions on Emerging Telecommunications, 2019;e3563. "[Code available here](https://github.com/victorcroisfelt/exp-lsf-spatial-corr-mmimo-chn-est)". 
- Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão, "Kaczmarz Precoding and Detection for Massive MIMO Systems", presented at IEEE Wireless Communications and Networking Conference (IEEE WCNC), Marrakech, Morocco on April 16, 2019. Publication is expected soon.

**About the code:** The codes provided herein can be used to simulate the Figs. 1 to 4 contained in the technical report (Appendix A of the FYP); this is done by running the scripts that have "simulation" in their names, while those with "function" in their names are called by the main scripts. Further details about each file can be found inside them.

## Citing this Repository and License
This code is licensed under the GPLv3 license. If you use any part of this repository for research, please consider to cite the FYP as:

Victor Croisfelt Rodrigues, "Spatial Correlation and Low Complexity Signal Processing Techniques in Massive MIMO Systems", Final Year Project, Universidade Estadual de Londrina, Londrina, Brazil, December, 2018.

