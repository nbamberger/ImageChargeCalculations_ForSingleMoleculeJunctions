function PlotAtomCoords(coordinates, atomic_numbers)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Plots a molecule by placing atomic symbols at
    %each of the atomic coordinates. 
    %
    %~~~INPUTS~~~:
    %
    %coordinates: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %atom_numbers: a vector listing the atomic number of each atom
    
    
    atomic_dictionary = importdata('atomic_symbol_dictionary.mat');
    atom_labels = atomic_dictionary(atomic_numbers);

    figure();
    hold on;
    for i = 1:size(coordinates,1)
        text(coordinates(i,1),coordinates(i,2),coordinates(i,3),atom_labels(i));
    end
    
    minL = min(min(coordinates))*1.1;
    maxL = max(max(coordinates))*1.1;
    xlim([minL maxL]);
    ylim([minL maxL]);
    zlim([minL maxL]);
    
    xlabel('X (nm)');
    ylabel('Y (nm)');
    zlabel('Z (nm)');

end