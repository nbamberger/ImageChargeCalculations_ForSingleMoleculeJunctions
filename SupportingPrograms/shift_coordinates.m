function new_coords = shift_coordinates(coords, shift_amount, axis)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: shifts a set of x,y,z coordinates by a given
    %amount along one of the cartesian axes (default is x)
    %
    %~~~INPUTS~~~:
    %
    %coords: a 3-column matrix of coordiantes with x, y, and z in the three
    %   columns, respectively, and one point in each row
    %
    %shift_amount: the amount to shift all coordinates by
    %
    %axis: the axis to rotate about; acceptable values are 'x', 'y', or
    %   'z'; 'x' is default if parameter is left out or empty
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %new_coords: shifted coordinates, again as a 3-column matrix
    
    
    if nargin < 3 || isempty(axis)
        axis = 'x';
    end
    
    %Make sure the specified axis makes sense:
    if ~strcmp(axis,'x') && ~strcmp(axis,'y') && ~strcmp(axis,'z')
        error("Unrecognized axis choice; use 'x', 'y', or 'z' ('x' is default if left blank or unspecified)");
    end
    
    if strcmp(axis,'z')
        dim = 3;
    elseif strcmp(axis,'y')
        dim = 2;
    elseif strcmp(axis,'x')
        dim = 1;
    end
    
    new_coords = coords;
    new_coords(:,dim) = new_coords(:,dim) + shift_amount;

end