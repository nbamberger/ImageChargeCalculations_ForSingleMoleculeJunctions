function [coordinates, interDist] = rotate_mol_to_align_axis(coordinates, ID1, ID2)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: rotate the coordinates of a molecule so that the
    %line between the atoms in positions ID1 and ID2 in the coordinate
    %matrix is aligned along the x-axis (and centered at the origin). Also
    %attempt to rotate the molecule to place it mostly in the xy-plane. 
    %
    %~~~INPUTS~~~:
    %
    %coordinates: 3-column array containing the x-, y-, and z-coordinates
    %   of each atom in a molecule
    %
    %ID1/ID2: the row IDs in coordinates that correspond to the two atoms
    %   that the user wants to align with the x-axis. Typically these would
    %   be the two binding atoms in a molecule (e.g. sulfurs, nitorgens).
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %coordinates: the new rotated (and shifted) coordinates
    %
    %interDist: the distance between the two atoms that have now been
    %   aligned along the x-axis
    

    current_vector = coordinates(ID2,:) - coordinates(ID1,:);

    %Rotate about z to align with x-axis
    ang = atan(current_vector(2)/current_vector(1))*180/pi;
    coordinates = rotate_coordinates(coordinates,-1*ang,'z');
    current_vector = coordinates(ID2,:) - coordinates(ID1,:);
    
    %Rotate about y to remove any tilt (seems this should NOT be negative;
    %always?  depends?  I'm not sure...)
    ang = atan(current_vector(3)/current_vector(1))*180/pi;
    coordinates = rotate_coordinates(coordinates,ang,'y');   

    %Remove any y- or z- offset
    coordinates(:,2) = coordinates(:,2) - coordinates(ID1,2);
    coordinates(:,3) = coordinates(:,3) - coordinates(ID1,3);
    
    %Shift so that molecule is centered along x between the two given atoms
    com = mean([coordinates(ID1,1),coordinates(ID2,1)]);
    coordinates(:,1) = coordinates(:,1) - com;
    
    %If molecules seems to be mostly in the xz-plane, rotate it 90 degrees
    %to put it mostly in the xy plane
    if std(coordinates(:,3)) > std(coordinates(:,2))
        coordinates = rotate_coordinates(coordinates, 90, 'x');
    end
    
    %Also calculate the distance between the two aligned atoms
    interDist = abs(coordinates(ID1,1) - coordinates(ID2,1));
  
end