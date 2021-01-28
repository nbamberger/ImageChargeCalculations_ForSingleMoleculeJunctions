function [GapSizes, Energies] = ...
    ImageChargeCalculation_GapSizeSweep_SphericalElectrodes(...
    charges, positions, radius, StartGap, EndGap, nSteps, minEnChange, maxIter, ToPlot)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Runs an image charge calculation for several 
    %gap sizes, assuming the two electrodes are each in the yz-plane and
    %symmetric about zero
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %positions: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %radius: the radius (in nm) of both spherical electrodes
    %
    %StartGap: the starting, smallest gap size, in units of nm.  Spherical
    %   electrodes will be centered on the x-axis, symmetrically around
    %   zero, such that the gap between them is this gap
    %
    %EndGap: the ending, largest gap size, in units of nm.  
    %
    %nSteps: the number of different gap sizes to run calculations for,
    %   which will be evenly spaced between StartGap and EndGap
    %
    %minEnChange: image charge iteration stops if the change in energy
    %   is less than this amount, in units of eV
    %
    %maxIter: image charge iteration will also stop if this number of
    %   iterations is reached
    %
    %ToPlot: logical variable; whether or not to make a plot showing how
    %   the renormalization energy varies with gap size
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %    
    %GapSizes: vector containing each of the gap sizes used
    %
    %Energies: vector containing the renormalization energy for each gap
    %   size
    
    
    %Default inputs
    if nargin < 9
        ToPlot = true;
    end
    if nargin < 8
        maxIter = 500;
    end
    if nargin < 7
        minEnChange = 0.005;
    end
    if nargin < 6
        nSteps = 10;
    end

    GapSizes = linspace(StartGap, EndGap, nSteps)';
    Energies = zeros(nSteps, 1);
    
    for i = 1:nSteps
        ens = SingleGapSize_RenormEnergy_SphericalElectrodes(charges,...
            positions,-GapSizes(i)/2-radius,GapSizes(i)/2+radius,radius,...
            radius,minEnChange,maxIter,false);
        
        %Average together last two energies since the energy oscillates as
        %it converges
        n = length(ens);
        Energies(i) = (ens(n-1) + ens(n))/2;
    end
    
    if ToPlot
        %Show system in the smallest gap size
        DisplaySystem_SphericalElectrodes(charges,positions,...
            -GapSizes(1)/2-radius,GapSizes(1)/2+radius,radius,radius,0);
        title('System at Smallest Gap Size');
        
        figure();
        plot(GapSizes, Energies,'-o');
        xlabel('Gap Size (nm)');
        ylabel('Image Charge Renormalization Energy (eV)');
    end

end