%This Matlab script can be used to generate Figure 1 in the technical 
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

%Number of UEs per BS
K = 18;

%% Scenario setup 

%Define BS radius (m)
cellRadius = 500;

%Distance between BSs (m)
interBSdistance = sqrt(3)*cellRadius;

%Define BS positions using complex coordinates (m)
BSlocations = [0+1i*0 interBSdistance*exp(1i*(pi/3*(0:5)))];
%Important: The center cell or home cell is considered to be the index 1,
%i.e., j = 1.

%Distribute UEs at random in the cell area
[UElocations] = functionDistributeUniformlyUEs(L,K,cellRadius,BSlocations);

%% Plotting the system setup

hFig = figure;
hold on

%Angle variation to plot a hexagon contour
theta = (pi/6):(pi/3):(14*pi/6);

%Go through all cells
for l = 1:L
    
    %Plotting the hexagon contour
    hexagonContour = cellRadius*(cos(theta)+1i*sin(theta)) + BSlocations(l);
    plot(real(hexagonContour),imag(hexagonContour),'k-','linewidth',2)
    
    %Plotting BS and UEs
    x = plot(real(BSlocations(l)),imag(BSlocations(l)),'ks','linewidth',2);
    y = plot(real(UElocations(:,l)),imag(UElocations(:,l)),'r*'); 
    
    %Check if the cell is the center one
    if l == 1
        
        %Shadding the center cell
        h =  fill(real(hexagonContour),imag(hexagonContour),'k','LineStyle','none');
        set(h,'facealpha',.1)
            
    end
    
end

%Legend of BS and UE
legend([x y],{'BS','UE'},'Location','best')

%Arranging some properties of the figure
color = get(hFig,'Color');

set(gca,'XColor',color,'YColor',color,'TickDir','out')

set(gca,'xtick',[])
set(gca,'ytick',[])