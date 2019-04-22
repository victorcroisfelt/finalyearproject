%This Matlab script can be used to generate Figure 4 in the technical
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

%Number of BSs
L = 7;

%Define the range of the number of BS antennas
Mrange = [10 100];

%Define the range of the number of UEs
Krange = 1:100;

%Extract the maximum values
Kmax = max(Krange);

%Number of setups with randomly generated statistics
nbrOfSetups = 1e3;

%% Scenario setup

%Define BS radius [m]
cellRadius = 500;

%Distance between BSs [m]
interBSdistance = sqrt(3)*cellRadius;

%Define BS positions using complex coordinates [m]
BSlocations = [0+1i*0 interBSdistance*exp(1i*(pi/3*(0:5)))]';

%% Propagation parameters

%Define total uplink transmit power per UE [mW]
rhoul = 200;

%Define total downlink transmit power per BS [mW]
rhodl = 200;

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

%Prepare to save the simulation results
sumSE_UL_MR = zeros(length(Krange),length(Mrange),nbrOfSetups);
sumSE_UL_ZF = zeros(length(Krange),length(Mrange),nbrOfSetups);
sumSE_DL_MR = zeros(length(Krange),length(Mrange),nbrOfSetups);
sumSE_DL_ZF = zeros(length(Krange),length(Mrange),nbrOfSetups);

%Go through all setups
for s = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(s) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Prepare to store pathloss numbers [dB]
    pathgaindB = zeros(Kmax,L);
    
    %Randomly distribute the UEs in the cell area
    UElocations = functionDistributeUniformlyUEs(L,Kmax,cellRadius,BSlocations);
    
    %Go through all cells
    for l = 1:L
        
        %Compute distances between UEs in BS l and BS j
        distancesBSj = abs(UElocations(:,l)-BSlocations(1));
        
        %Compute distant-dependent path gains [dB]
        pathgaindB(:,l) = -alpha*10*log10(distancesBSj);
        
    end
    
    %Compute the normalized channel gains, where the normalization is by
    %the noise power
    channelGaindB = pathgaindB-noiseVariancedBm;
    
    %Compute the large-scale coefficients, adding the shadowing term
    betasLoop = 10.^(std_shad.*randn(Kmax,L)./10).*10.^(channelGaindB./10);
    
    %Go through all number of UEs
    for k = 1:length(Krange)
        
        %Extract the current values of some variables
        K = Krange(k);
        betas = betasLoop(1:K,:);
        
        %Go through all number of BS antennas
        for m = 1:length(Mrange)
            
            %Extract the current values of some variables
            M = Mrange(m);
            
            %MMSE channel estimation
            [taup,psis] = functionChannelEstimates(L,K,rhoul,betas);
            
            %Compute UL SE
            [SE_UL_MR,SE_UL_ZF] = functionComputeSE_UL(L,M,K,tauc,taup,rhoul,betas,psis);
            
            %Compute DL SE
            [SE_DL_MR,SE_DL_ZF] = functionComputeSE_DL(L,M,K,tauc,taup,rhodl,betas,psis);
            
            %Store the simulation results
            sumSE_UL_MR(k,m,s) = sum(SE_UL_MR);
            sumSE_UL_ZF(k,m,s) = sum(SE_UL_ZF);
            sumSE_DL_MR(k,m,s) = sum(SE_DL_MR);
            sumSE_DL_ZF(k,m,s) = sum(SE_DL_ZF);
            
        end
        
    end
    
end

%% Plot the simulation results

%UL: Average sum SE vs Number of UEs
figure;
subplot(1,2,1)
hold on; box on; ax = gca;

%Compute the average result for MR
y10_MR = mean(sumSE_UL_MR(:,1,:),3); % M = 10
y100_MR = mean(sumSE_UL_MR(:,2,:),3); % M = 100

%Compute the average result for ZF
y10_ZF = mean(sumSE_UL_ZF(:,1,:),3); % M = 10
y100_ZF = mean(sumSE_UL_ZF(:,2,:),3); % M = 100

%Plot referece markers
plot(Krange(1),y10_MR(1),'d--','LineWidth',1)
plot(Krange(1),y10_ZF(1),'--','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot referece markers
plot(Krange(1),y100_MR(1),'d-.','LineWidth',1)
plot(Krange(1),y100_ZF(1),'-.','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot the lines
plot(Krange(1:Mrange(1)),y10_MR(1:Mrange(1)),'--','LineWidth',1)
plot(Krange(1:Mrange(1)),y10_ZF(1:Mrange(1)),'--','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot the lines
plot(Krange(1:Mrange(2)),y100_MR(1:Mrange(2)),'--','LineWidth',1)
plot(Krange(1:Mrange(2)),y100_ZF(1:Mrange(2)),'--','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot the missing markers
plot(Krange(5:5:Mrange(1)),y10_MR(5:5:Mrange(1)),'d','LineWidth',1)
ax.ColorOrderIndex = 1;
plot(Krange(10:10:Mrange(2)),y100_MR(10:10:Mrange(2)),'d','LineWidth',1)

xlabel('Number of UEs ($K$)');
ylabel('Average sum UL SE [bit/s/Hz]');

legend('MR, $M = 10$','ZF, $M = 10$','MR, $M = 100$','ZF, $M = 100$','Location','best');

%DL: Average sum SE vs Number of UEs
subplot(1,2,2)
hold on; box on; ax = gca;

%Compute the average result for MR
y10_MR = mean(sumSE_DL_MR(:,1,:),3); % M = 10
y100_MR = mean(sumSE_DL_MR(:,2,:),3); % M = 100

%Compute the average result for ZF
y10_ZF = mean(sumSE_DL_ZF(:,1,:),3); % M = 10
y100_ZF = mean(sumSE_DL_ZF(:,2,:),3); % M = 100

%Plot referece markers
plot(Krange(1),y10_MR(1),'d--','LineWidth',1)
plot(Krange(1),y10_ZF(1),'--','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot referece markers
plot(Krange(1),y100_MR(1),'d-.','LineWidth',1)
plot(Krange(1),y100_ZF(1),'-.','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot the lines
plot(Krange(1:Mrange(1)),y10_MR(1:Mrange(1)),'--','LineWidth',1)
plot(Krange(1:Mrange(1)),y10_ZF(1:Mrange(1)),'--','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot the lines
plot(Krange(1:Mrange(2)),y100_MR(1:Mrange(2)),'--','LineWidth',1)
plot(Krange(1:Mrange(2)),y100_ZF(1:Mrange(2)),'--','LineWidth',1)

%Reset the colors
ax.ColorOrderIndex = 1;

%Plot the missing markers
plot(Krange(5:5:Mrange(1)),y10_MR(5:5:Mrange(1)),'d','LineWidth',1)
ax.ColorOrderIndex = 1;
plot(Krange(10:10:Mrange(2)),y100_MR(10:10:Mrange(2)),'d','LineWidth',1)

xlabel('Number of UEs ($K$)');
ylabel('Average sum UL SE [bit/s/Hz]');

legend('MR, $M = 10$','ZF, $M = 10$','MR, $M = 100$','ZF, $M = 100$','Location','best');