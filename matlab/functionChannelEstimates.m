function [taup,psis,psis_singlecell,NMSE,NMSE_singlecell] = functionChannelEstimates(L,K,rhoul,betas)
%Compute the normalized mean square error (NMSE) and the average estimated
%variance for a canonical M-MIMO system operating in TDD mode. Also, those
%parameters are computed for a single cell case, where only the cell of
%interest (cell j) is considered to be operating. All the metrics are in
%relation to BS j.
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
%   K: number of UEs inside eachs cell.
%   rhoul: Uplink transmit power per UE (same for everyone).
%   betas: K x L matrix with the average large-scale coefficient of the
%   users over the entire system in relation to the center cell (cell j).
%
%@Outputs:
%   taup: Pilot length (samples).
%   psis: K x L matrix with the average variance of the channel estimates
%   over the entire system in relation to the center cell.
%   psis_singlecell: K x 1 vector with the average variance of the channel
%   estimates over the entire system in relation to the center cell for the
%   case considering only the operation of cell j.
%   NMSE: K x L matrix with the NMSE of each user in relation to cell j.
%   NMSE_singlecell: K x 1 vector with the NMSE of each user in relation to
%   cell j for the single-cell case.
%
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

%Length of pilot sequences [samples]
taup = K;

%Compute the normalized pilot transmitted power (W)
rhop = taup*rhoul;

%Prepare to store the simulation results
psis = zeros(K,L);

%Check the number of outputs
if nargout > 2
    
    psis_singlecell = zeros(K,1);
    
end

%Check the number of outputs
if nargout > 3
    
    NMSE = zeros(K,L);
    NMSE_singlecell = zeros(K,1);
    
end

%Go through all UEs in cell j
for k = 1:K
    
    %Compute the matrix denominator
    denom = 1 + rhop*sum(betas(k,:));
    
    if nargout > 2
        
        %Compute the matrix denominator
        denom_singlecell = 1+rhop*betas(k,1);
        
        %Compute the average estimated variance for the single-cell case
        psis_singlecell(k) = (rhop*betas(k,1)^2)/denom_singlecell;
        
    end
    
    %Go through all cells
    for l = 1:L
        
        %Compute the average estimated variance
        psis(k,l) = (rhop*betas(k,l)^2)/denom;
        
        if nargout >= 3
            
            %Compute the NMSE per user per antenna
            NMSE(k,l) = (betas(k,l)-psis(k,l))/betas(k,l);
            
        end
        
    end
    
    if nargout >= 3
        
        %Compute the NMSE per user per antenna for the single-cell case
        NMSE_singlecell(k) = (betas(k,1)-psis_singlecell(k))/betas(k,1);
        
    end
    
end
