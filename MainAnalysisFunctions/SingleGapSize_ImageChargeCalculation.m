function energy = SingleGapSize_ImageChargeCalculation(charges, positions, ...
    left_plane, right_plane, minEnChange, maxIter, ToPlot)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Calculates the image charge renormalization 
    %energy for a single position of a set of partial charges between two
    %planar electrodes in the yz plane.  Returns a vector of energies after
    %each iteration. Keeps iterating until the change in energy is less 
    %than minEnChange.
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %positions: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %left_plane/right_plane: the y-coordinate of the left/right planar
    %   electrode, respectively (assumed to be in the yz-plane), in units
    %   of nm
    %
    %minEnChange: image charge iteration stops if the change in energy
    %   is less than this amount, in units of eV
    %
    %maxIter: image charge iteration will also stop if this number of
    %   iterations is reached
    %
    %ToPlot: logical variable; whether or not to make a plot showing how
    %   the cumulative energy changed after each iteration
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %    
    %energy: a vector containing the cumulative image charge interaction
    %   energy after each iteration
    

    %Default inputs
    if nargin < 5
        minEnChange = 0.005;
    end
    if nargin < 6
        maxIter = 500;
    end
    if nargin < 7
        ToPlot = false;
    end

    energy = zeros(maxIter, 1);
    
    %Display a warning if any of the partial charges are on the wrong side
    %of the image planes
    if min(positions(:,1)) < left_plane || max(positions(:,1)) > right_plane
        warning('Atom(s) are inside electrodes!!');
    end
    
    %Calculate the energy after the first iteration
    energy(1) = interaction_energy_with_jth_images(charges, positions, ...
        left_plane, right_plane, 1);
    
    %Keep iterating until energy change is less than specified value or
    %maximum # of iterations is reached
    counter = 1;
    EnChange = Inf;
    while abs(EnChange) > minEnChange && counter < maxIter
        counter = counter + 1;
        EnChange = interaction_energy_with_jth_images(charges, positions, ...
            left_plane, right_plane, counter);
        energy(counter) = energy(counter - 1) + EnChange;
    end
    energy = energy(1:counter);
    
    %Display warning if max iterations were reached:   
    if abs(EnChange) > minEnChange
        warning('Reached maximum # of iterations without meeting convergence criteria');
        disp(strcat('Threshhold:', {' '}, num2str(minEnChange), ' Last Change:', ...
            {' '}, num2str(EnChange)));
    end
    
    if ToPlot
        energyPlot = [0; energy];
        figure();
        plot((0:counter)',energyPlot);
        xlabel('Iteration #');
        ylabel('Interaction Energy (eV)');
    end

end