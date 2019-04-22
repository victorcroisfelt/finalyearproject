function [SE_UL_MR,SE_UL_ZF,SE_UL_MR_singlecell,SE_UL_ZF_singlecell] = functionComputeSE_UL(L,M,K,tauc,taup,rhoul,betas,psis,psis_singlecell)
%Evaluate the uplink spectral efficiency (SE) through the application of
%the maximum-ratio (MR) and zero-forcing (ZF) combining schemes, when
%considering a canonical multicell Massive MIMO system operating in TDD.
%Also, the single-cell case is computed as a baseline.
%
%This Matlab function is used in the technical report - "Massive MIMO
%System in TDD Mode: Channel Estimation and Spectral Efficiency" - included
%in the following final year project (FYP):
%
%Victor Croisfelt Rodrigues, "Spatial Correlation and Low Complexity Signal
%Processing Techniques in Massive MIMO Systems", Final Year Project,
%Universidade Estadual de Londrina, Londrina, Brazil, December, 2018.
%
%Download FYP: https://github.com/victorcroisfelt/finalyearproject
%
%This is version 3.0 (Last edited: 04-21-2019)
%
%License: This code is licensed under the GPLv3 license. If you in any way
%use this code for research that results in publications, please reference
%our original FYP as shown above.
%
%@Inputs:
%   L: Number of BSs or cells.
%   M: Number of antennas per BS.
%   K: number of UEs inside eachs cell.
%   tauc: Length of the coherence block.
%   taup: Pilot length (samples).
%   rhoul: Uplink transmit power per UE (same for everyone).
%   betas: K x L matrix with the average large-scale coefficient of the
%   users over the entire system in relation to the center cell (cell j).
%   psis: K x L matrix with the average variance of the channel estimates
%   over the entire system in relation to the center cell.
%   psis_singlecell: K x 1 vector with the average variance of the channel
%   estimates over the entire system in relation to the center cell for the
%   case considering only the operation of cell j.
%
%@Outputs:
%   SE_UL_MR: K x 1 vector where the kth element is the uplink SE of UE k
%   in cell j achieved with MR combining.
%   SE_UL_ZF: K x 1 vector where the kth element is the uplink SE of UE k
%   in cell j achieved with ZF combining.
%   SE_UL_MR_singlecell: K x 1 vector where the kth element is the uplink
%   SE of UE k in cell j achieved with MR combining for the single-cell
%   case.
%   SE_UL_ZF_singlecell: K x 1 vector where the kth element is the uplink
%   SE of UE k in cell j achieved with ZF combining for the single-cell
%   case.
%
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

%Compute the pre-log factor
prelogFactor = (0.5)*(tauc-taup)/(tauc);

%Prepare to store simulation results for the multi-cell case
SE_UL_MR = zeros(K,1);
SE_UL_ZF = zeros(K,1);

%Prepare to store the simulation results for single-cell case
if nargin > 8
    
    SE_UL_MR_singlecell = zeros(K,1);
    SE_UL_ZF_singlecell = zeros(K,1);
    
end

%Go through all UEs in cell j
for k = 1:K
    
    %Multicell
    
    %Maximum-Ratio (MR)
    signal = M*rhoul*psis(k,1);
    interfnoise = 1+rhoul*sum(betas(k,:))+M*rhoul*sum(psis(k,2:L));
    
    SE_UL_MR(k) = prelogFactor*real(log2(1+signal/interfnoise));
    
    %Zero-Forcing (ZF)
    signal = (M-K)*rhoul*psis(k,1);
    interfnoise = 1+rhoul*sum(betas(k,:)-psis(k,:))+(M-K)*rhoul*sum(psis(k,2:L));
    
    SE_UL_ZF(k) = prelogFactor*real(log2(1+signal/interfnoise));
    
    %Single-cell
    if nargin > 8
        
        %Maximum-Ratio (MR)
        signal = M*rhoul*psis_singlecell(k);
        interfnoise = 1+rhoul*sum(betas(k,1));
        
        SE_UL_MR_singlecell(k) = prelogFactor*real(log2(1+signal/interfnoise));
        
        %Zero-Forcing (ZF)
        signal = (M-K)*rhoul*psis_singlecell(k);
        interfnoise = 1+rhoul*sum(betas(k,1)-psis_singlecell(k));
        
        SE_UL_ZF_singlecell(k) = prelogFactor*real(log2(1+signal/interfnoise));
        
    end
    
end