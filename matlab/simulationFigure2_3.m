%This Matlab script can be used to generate Figs. 2 and 3 in the technical
%report - "Massive MIMO System in TDD Mode: Channel Estimation and Spectral
%Efficiency" - included in the following final year project (FYP):
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
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

%Initialization
close all;
clearvars;

%% Simulation parameters

%Choose the desired simulation subfigure:
%   simulation == 1: Figure 2
%   simulation == 2: Figure 3
%
simulation = 1;

%Number of BSs
L = 7;

%Define parameters in each simulation scenario
if simulation == 1
    
    %Number of BS antennas
    M = 100;
    
    %Determine maximum number of BS antennas
    Mmax = M;
    
    %Number of UEs per BS
    K = 18;
    
    %Select range of total uplink transmit power per UE [mW]
    rhoul_range = logspace(-3,2,10);
    
    %Number of measurement points along horizontal axis in final figure
    nbrOfPoints = length(rhoul_range);
    
    %Number of setups with randomly generated statistics
    nbrOfSetups = 1e3;
    
elseif simulation == 2
    
    %Select range of number of BS antennas
    Mrange = round(logspace(2,5,10));
    
    %Determine maximum number of BS antennas
    Mmax = max(Mrange);
    
    %Number of UEs per BS
    K = 18;
    
    %Define total uplink transmit power per UE [mW]
    rhoul = 200;
    
    %Define total downlink transmit power per BS [mW]
    rhodl = 200;
    
    %Number of measurement points along horizontal axis in final figure
    nbrOfPoints = length(Mrange);
    
    %Number of setups with randomly generated statistics
    nbrOfSetups = 1e3;
    
end

%% Scenario setup

%Define BS radius [m]
cellRadius = 500;

%Distance between BSs [m]
interBSdistance = sqrt(3)*cellRadius;

%Define BS positions using complex coordinates [m]
BSlocations = [0+1i*0 interBSdistance*exp(1i*(pi/3*(0:5)))]';
%Important: The center cell or home cell is considered to be the index 1,
%i.e., j = 1.

%% Propagation parameters

%Pathloss exponent
alpha = 3.76;

%Standard deviation of large-scale fading variations (shadowing) [dB]
std_shad = 8;

%Communication bandwidth [Hz]
bandwidth = 20e6;

%Define noise figure at BS [dB]
noiseFigure = 10;

%Compute total noise power [dBm]
noiseVariancedBm = -174+10*log10(bandwidth)+noiseFigure;

%Select length of coherence block [samples]
tauc = 200;

%% Simulation

%Prepare to save simulation results
if simulation == 1
    
    %Channel estimation error based on NMSE (normalized mean squared error)
    meanNMSEj = zeros(nbrOfPoints,nbrOfSetups); % within cell j
    meanNMSEl = zeros(nbrOfPoints,nbrOfSetups); % surrounding cells
    meanNMSE_singlecell = zeros(nbrOfPoints,nbrOfSetups); % only cell j
    
elseif simulation == 2
    
    %Uplink (UL) spectral efficiency (SE) per user
    meanSE_UL_MR = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_UL_ZF = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_UL_MR_singlecell = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_UL_ZF_singlecell = zeros(nbrOfPoints,nbrOfSetups);
    
    %Downlink (DL) spectral efficiency (SE) per user
    meanSE_DL_MR = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_DL_ZF = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_DL_MR_singlecell = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_DL_ZF_singlecell = zeros(nbrOfPoints,nbrOfSetups);
    
end

%Go through all setups
for s = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(s) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Prepare to store pathloss numbers [dB]
    pathgaindB = zeros(K,L);
    
    %Randomly distribute the UEs in the cell area
    UElocations = functionDistributeUniformlyUEs(L,K,cellRadius,BSlocations);
    
    %Go through all cells
    for l = 1:L
        
        %Compute distances between UEs in BS l and BS j
        distancesBSj = abs(UElocations(:,l)-BSlocations(1));
        
        %Compute distant-dependent path gains [dB]
        pathgaindB(:,l) = -alpha*10*log10(distancesBSj);
        
    end
    
    %Compute the normalized channel gains, where the normalization is by
    %the noise power
    channelGaindB = pathgaindB - noiseVariancedBm;
    
    %Compute the large-scale coefficients, adding the shadowing term
    betas = 10.^(std_shad.*randn(K,L)./10).*10.^(channelGaindB./10);
    
    %Go through all points
    for m = 1:nbrOfPoints
        
        %Output simulation progress
        disp([num2str(m) ' points out of ' num2str(nbrOfPoints)]);
        
        %Check the simulation scenario
        if simulation == 1
            
            %Extract current uplink transmitted power
            rhoul = rhoul_range(m);
            
        elseif simulation == 2
            
            %Extract current number of BS antennas
            M = Mrange(m);
            
        end
        
        if simulation == 1
            
            %MMSE channel estimation
            [taup,psis,psis_singlecell,NMSE,NMSE_singlecell] = functionChannelEstimates(L,K,rhoul,betas);
            
            %Store the simulation results
            meanNMSEj(m,s) = mean(NMSE(:,1));
            meanNMSEl(m,s) = mean(mean(NMSE(:,2:L),1),2);
            meanNMSE_singlecell(m,s) = mean(NMSE_singlecell);
            
        elseif simulation == 2
            
            %MMSE channel estimation
            [taup,psis,psis_singlecell] = functionChannelEstimates(L,K,rhoul,betas);
            
            %Compute UL SE
            [SE_UL_MR,SE_UL_ZF,SE_UL_MR_singlecell,SE_UL_ZF_singlecell] = functionComputeSE_UL(L,M,K,tauc,taup,rhoul,betas,psis,psis_singlecell);
            
            %Compute DL SE
            [SE_DL_MR,SE_DL_ZF,SE_DL_MR_singlecell,SE_DL_ZF_singlecell] = functionComputeSE_DL(L,M,K,tauc,taup,rhodl,betas,psis,psis_singlecell);
            
            %Store the simulation results
            
            %UL
            meanSE_UL_MR(m,s) = mean(SE_UL_MR);
            meanSE_UL_ZF(m,s) = mean(SE_UL_ZF);
            meanSE_UL_MR_singlecell(m,s) = mean(SE_UL_MR_singlecell);
            meanSE_UL_ZF_singlecell(m,s) = mean(SE_UL_ZF_singlecell);
            
            %DL
            meanSE_DL_MR(m,s) = mean(SE_DL_MR);
            meanSE_DL_ZF(m,s) = mean(SE_DL_ZF);
            meanSE_DL_MR_singlecell(m,s) = mean(SE_DL_MR_singlecell);
            meanSE_DL_ZF_singlecell(m,s) = mean(SE_DL_ZF_singlecell);
            
        end
        
    end
    
end

%% Plot simulation results

if simulation == 1
    
    figure;
    hold on; box on;
    
    plot(rhoul_range,10*log10(mean(meanNMSE_singlecell,2)),'k:','LineWidth',1.25)
    plot(rhoul_range,10*log10(mean(meanNMSEj,2)),'-','LineWidth',1)
    plot(rhoul_range,10*log10(mean(meanNMSEl,2)),'o--','LineWidth',1)
    
    ylabel('Average NMSE [dB]');
    xlabel('Radiated power per terminal [mW]');
    
    legend('Single-cell','Multi-cell: home cell','Multi-cell: other cells','Location','SouthWest');
    set(gca,'XScale','log');
    
elseif simulation == 2
    
    figure;
    
    %Uplink
    subplot(1,2,1)
    hold on; box on;
    
    plot(Mrange,(mean(meanSE_UL_MR_singlecell,2)),'k:s','LineWidth',1.25)
    plot(Mrange,(mean(meanSE_UL_ZF_singlecell,2)),'k:h','LineWidth',1.25)
    plot(Mrange,(mean(meanSE_UL_MR,2)),'--s','LineWidth',1)
    plot(Mrange,(mean(meanSE_UL_ZF,2)),'--h','LineWidth',1)
    
    ylabel('Average UL SE [bit/s/Hz/user]');
    xlabel('Number of antennas ($M$)');
    
    legend('Single-cell: MR','Single-cell: ZF','Multi-cell: MR','Multi-cell: ZF','Location','NorthWest');
    
    set(gca,'XScale','log');
    
    ylim([0 11])
    
    %Downlink
    subplot(1,2,2)
    hold on; box on;
    
    plot(Mrange,(mean(meanSE_DL_MR_singlecell,2)),'k:s','LineWidth',1.25)
    plot(Mrange,(mean(meanSE_DL_ZF_singlecell,2)),'k:h','LineWidth',1.25)
    plot(Mrange,(mean(meanSE_DL_MR,2)),'--s','LineWidth',1)
    plot(Mrange,(mean(meanSE_DL_ZF,2)),'--h','LineWidth',1)
    
    ylabel('Average DL SE [bit/s/Hz/user]');
    xlabel('Number of antennas ($M$)');
    
    legend('Single-cell: MR','Single-cell: ZF','Multi-cell: MR','Multi-cell: ZF','Location','NorthWest');
    
    set(gca,'XScale','log');
    
    ylim([0 11])
    
end
