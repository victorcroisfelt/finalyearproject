function [UElocations] = functionDistributeUniformlyUEs(L,K,cellRadius,BSlocations)
%Uniformly drops the UEs within a hexagon cell shape. The hexagon cell has
%an radius equals to cellRadius with the center located at the BSlocations. 
%All the distances are computed in meters [m].
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
%   K: number of UEs inside the cell.
%   cellRadius: cell radius (m).
%   BSlocations: L x 1 vector with the BS position.
%
%@Outputs:
%   UElocations : K x L matrix with UEs positions in meters at each BS.
%
%References:
%[1] Marzetta, T., Larsson, E., Yang, H., & Ngo, H. (2016). "Fundamentals 
%of Massive MIMO". Cambridge: Cambridge University Press. 
%doi:10.1017/CBO9781316799895
%

    %Prepare to store the locations of each user
    UElocations = zeros(K,L); 

    %Go through all cells
    for ll = 1:L
        
        % Go through all UEs
        for kk = 1:K
            
            %Drawing the sector
            sector = rand;

            %Square uniform distribution
            u = cellRadius*rand;
            v = cellRadius*rand;

            %Mapping u and v conform the sector drawn
            if (sector > 0 && sector < 1/3)
                
                UElocations(kk,ll) = (sqrt(3)/2*u) + 1i*(-1/2*u+v) + BSlocations(ll); 

            elseif (sector >= 1/3 && sector < 2/3)
                
                UElocations(kk,ll) = (-sqrt(3)/2*u + sqrt(3)/2*v) + 1i*(-1/2*u-1/2*v) + BSlocations(ll);

            else
                
                UElocations(kk,ll) = (-sqrt(3)/2*v) + 1i*(u-1/2*v) + BSlocations(ll);

            end
            
        end
        
    end
    
end