function Energies = SingleGapSize_RenormEnergy_SphericalElectrodes(...
    charges, positions, left_center, right_center, left_radius, ...
    right_radius, minEnchange, maxIter,ToPlot)
    %Copyright 2021 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Iteratively calculates the image charge
    %interaction energy for a set of charges placed between two spherical
    %electrodes, recording the cumulative energy after each successive set
    %of image charges is added and proceeding until convergence criteria is
    %reached.
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %positions: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %left_center/right_center: the x-coordinate center (in nm) of the 
    %   left/right spherical electrode; the y- and z-coordinates of both
    %   centers are assumed to be zero, i.e. both centers lie on the x-axis
    %
    %left_radius/right_radius: the radius (in nm) of the left/right
    %   spherical electrode
    %
    %minEnChange: image charge iteration stops if the change in energy
    %   is less than this amount, in units of eV
    %
    %maxIter: image charge iteration will also stop if this number of
    %   iterations is reached
    %
    %ToPlot: logical variable; whether or not to make a plot showing how
    %   the renormalization energy varies as more and more image charges
    %   are added
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %
    %Energies: vector containing the renormalization energy after each
    %   successive set of image charges is added; average together the last
    %   two values to get the best approximation for the true
    %   renormalization energy
    

    if nargin < 9
        ToPlot = false;
    end
    
    %Minimum # of iterations
    minIter = 5;

    %Define some constants:
    coulomb_constant = 8987551787.3681764; %in N m^2/C^2
    elementary_charge = 1.602176634E-19; %in C
    joules_to_eV = 6.242E+18;
    TotalConstant = coulomb_constant*elementary_charge^2*1E9*joules_to_eV;
    
    Energies = zeros(maxIter,1);
    n = length(charges);

    %We will explicitly compute the first iteration; get the charges and
    %positions of this first set of images
    [left_ch, left_co] = getSphere_ImageCharges(charges, positions, ...
        left_center, left_radius);
    [right_ch, right_co] = getSphere_ImageCharges(charges, positions, ...
        right_center, right_radius);    
    
    %Now add up all the energies for the original charges interacting with
    %the first set of image charges
    E = 0;
    for i = 1:n
        for j = 1:n
            E = E + charges(i)*left_ch(j)/(2*sqrt((positions(i,1)-left_co(j,1))^2 + ...
                (positions(i,2)-left_co(j,2))^2 + (positions(i,3)-left_co(j,3))^2)) + ...
                charges(i)*right_ch(j)/(2*sqrt((positions(i,1)-right_co(j,1))^2 + ...
                (positions(i,2)-right_co(j,2))^2 + (positions(i,3)-right_co(j,3))^2));
        end
    end
    Energies(1) = E*TotalConstant;
    
    %Now proceed with the 2nd, etc., iterations
    iter = 2;
    enChange = E;
    while (abs(enChange) > minEnchange && iter <= maxIter) || iter < minIter
        
        %Get the charges and positions of the next order of image charges
        left_ch_old = left_ch;
        left_co_old = left_co;      
        [left_ch, left_co] = getSphere_ImageCharges(right_ch, right_co, ...
            left_center, left_radius);
        [right_ch, right_co] = getSphere_ImageCharges(left_ch_old, left_co_old, ...
            right_center, right_radius); 
        
        %Add up all the energies for the original charges interacting with
        %this next set of image charges
        E = 0;
        for i = 1:n
            for j = 1:n
                E = E + charges(i)*left_ch(j)/(2*sqrt((positions(i,1)-left_co(j,1))^2 + ...
                    (positions(i,2)-left_co(j,2))^2 + (positions(i,3)-left_co(j,3))^2)) + ...
                    charges(i)*right_ch(j)/(2*sqrt((positions(i,1)-right_co(j,1))^2 + ...
                    (positions(i,2)-right_co(j,2))^2 + (positions(i,3)-right_co(j,3))^2));
            end
        end  
        
        %Store cumulative energy after this iteration
        enChange = E*TotalConstant;
        Energies(iter) = Energies(iter-1) + enChange;
        iter = iter + 1;
        
    end
    
    %Display warning if convergence not achieved; otherwise trim energy
    %vector
    if iter == maxIter + 1
        warning('Maximum # of iterations reached!');
    else
        Energies = Energies(1:iter-1);
    end
    
    %Make plot if requested
    if ToPlot
        energyPlot = [0; Energies];
        figure();
        plot((0:iter-1)',energyPlot,'-o');
        xlabel('Iteration #');
        ylabel('Cumulative Interaction Energy (eV)');
    end

end