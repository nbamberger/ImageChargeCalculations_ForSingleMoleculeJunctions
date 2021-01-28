function [sphere_ch, sphere_co] = getSphere_ImageCharges(charges, coords, ...
    sphere_center, sphere_radius)
    %Copyright 2021 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Given a set of charges and coordinates, 
    %calculate the new charges and coordinates for the images of those
    %point charges inside a conductive sphere. The sphere is assumed to be
    %centered at a point on the x-axis, so only a single value (the 
    %x-coordinate) is needed to specify its center.
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %coords: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %sphere_center: the x-coordinate (in nm) of the center of the spherical
    %   electrode; y- and z-coordinates assumed to be zero
    %
    %sphere_radius: the radius (in nm) of the spherical electrode
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %    
    %sphere_ch: vector of charges (in units of e) for each image charge
    %   inside the spherical electrode
    %
    %sphere_co: 3-column matrix of the x-, y-, and z-coordinates (in nm)
    %   for each image charge inside the spherical electrode
    

    %Temporarily move coordinate system to put sphere center at center
    coords(:,1) = coords(:,1) - sphere_center;
    
    %Calculate distances of each point charge from sphere center
    dists = sqrt(coords(:,1).^2 + coords(:,2).^2 + coords(:,3).^2);
    
    %Display warning if any atoms fall inside the electrodes
    if min(dists) < sphere_radius
        warning('One or more atoms is inside the electrodes!');
    end
    
    %Find new charge amounts
    sphere_ch = -charges * sphere_radius ./ dists;
    
    %Find what new distances from sphere center should be, use to set new
    %coordinates
    new_dists = sphere_radius^2 ./ dists;   
    sphere_co = coords;
    for i = 1:3
        sphere_co(:,i) = sphere_co(:,i) .* new_dists ./ dists;
    end
    
    %Move coordinate system back to where it was
    sphere_co(:,1) = sphere_co(:,1) + sphere_center;

end