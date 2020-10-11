function energies = ImageRenorm_DuringRotation(charges, positions, ...
    initial_gap, angle_list, changeGapSize, ToPlot)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: rotates a given molecule about the z-axis to
    %multiple different rotation angles, and calculates the image charge
    %interaction energy at each rotation angle.  
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %positions: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %initial_gap: the (initial) gap between the two electrodes (will be 
    %   modified if changeGapSize is true).  Both electrodes are assumed to
    %   be in the yz-plane and symmetric about zero.  So, for example, a
    %   gap size of 1 means that the left electrode is at x = -0.5 nm and
    %   the right electrode is at x = 0.5 nm.  
    %
    %angle_list: a vector containing each angle to consider, in degrees
    %
    %changeGapSize: logical variable; if true, then the "initial_gap" is
    %   assumed to be the gap size at zero rotation angle, and for every
    %   other angle the gap size is adjusted to keep the spacing between
    %   molecule and electrodes roughly the same.  
    %
    %ToPlot: logical variable; whether or not to create a plot showing
    %   image charge renormalization energy vs. rotation angle
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %  
    %energies: vector containing the renormalization energy for each
    %   rotation angle considered
    
    
    %Default inputs
    if nargin < 5
        changeGapSize = false;
    end
    if nargin < 6
        ToPlot = true;
    end

    Nangles = length(angle_list);
    energies = zeros(Nangles, 1);
    gap_sizes = zeros(Nangles, 1);
    
    for i = 1:Nangles
        disp([i Nangles]);
        
        ang = angle_list(i) * pi/180;
        
        %Figure out gap size
        if changeGapSize
            gap_sizes(i) = initial_gap * cos(ang);
        else
            gap_sizes(i) = initial_gap;
        end
        
        newcoords = rotate_coordinates(positions,angle_list(i),'z');
        
        [~, energies(i)] = ImageChargeCalculation_GapSizeSweep(charges, ...
            newcoords, gap_sizes(i), gap_sizes(i), 1, 0.005, 500, false);
        
    end
    
    if ToPlot
        figure();
        plot(angle_list, energies, '-o');
        xlabel('Angle (degrees)');
        ylabel('Energy (eV)');
    end

end