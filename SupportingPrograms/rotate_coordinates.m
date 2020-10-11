function new_coords = rotate_coordinates(coords, rot_angle_degs, axis)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: rotates a set of x,y,z coordinates by a given
    %angle about one of the cartesian axes (default is z)
    %
    %~~~INPUTS~~~:
    %
    %coords: a 3-column matrix of coordiantes with x, y, and z in the three
    %   columns, respectively, and one point in each row
    %
    %rot_angle_degs: angle to rotate by, in degrees
    %
    %axis: the axis to rotate about; acceptable values are 'x', 'y', or
    %   'z'; 'z' is default if parameter is left out or empty
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %new_coords: rotated coordinates, again as a 3-column matrix
    
    
    if nargin < 3 || isempty(axis)
        axis = 'z';
    end
    
    %Make sure the specified axis makes sense:
    if ~strcmp(axis,'x') && ~strcmp(axis,'y') && ~strcmp(axis,'z')
        error("Unrecognized axis choice; use 'x', 'y', or 'z' ('z' is default if left blank or unspecified)");
    end

    %Convert angle to radians
    ang = rot_angle_degs * pi()/180;
    
    if strcmp(axis,'z')
        rot_matrix = [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1];
    elseif strcmp(axis,'y')
        rot_matrix = [cos(ang) 0 sin(ang); 0 1 0; -sin(ang) 0 cos(ang)];
    elseif strcmp(axis,'x')
        rot_matrix = [1 0 0; 0 cos(ang) -sin(ang); 0 sin(ang) cos(ang)];
    end
    
    %Perform the rotation (need to mess with the orientation of the
    %matrices first and after)
    coords = coords';
    new_coords = rot_matrix*coords;
    new_coords = new_coords';

end