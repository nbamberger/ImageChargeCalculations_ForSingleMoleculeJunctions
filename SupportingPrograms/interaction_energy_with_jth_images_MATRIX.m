function Emat = interaction_energy_with_jth_images_MATRIX(charges, positions, ...
    left_plane, right_plane, iterNum)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Given a set of partial charges (in e) and their
    %positions (in nm) and the positions of two image planes (assumed to be
    %in the yz-plane, and so specified just by their x-coordinates), this
    %function calculates the interaction energy between the original 
    %partial charges and their jth images on both the left and the right.
    %Returns this interaction energy as a matrix so that you can see how 
    %much each pair of partial charges contributes.
    %
    %~~~INPUTS~~~:
    %
    %charges: vector listing partial charges in units of e
    %
    %positions: 3-column array listing x-, y-, and x-coordinates of each
    %   partial charge in units of nm
    %
    %left_plane/right_plane: x-coordinate of left/right image plane in nm
    %
    %iterNum: what order image charges are being considered, A.K.A j of
    %   "jth images"
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %Emat: matrix in which the (i,k)th element represents the interaction
    %   energy between the ith original partial charge and the jth images
    %   of the kth partial charge in both the left and right electrodes, in
    %   units of eV
    

    %Find the coordinates of the images on the left (only need to change
    %x-coordinates);
    left_coords = positions;
    left_coords(:,1) = jth_left_image_Xcoords(left_coords(:,1), left_plane, ...
        right_plane, iterNum);
    
    %Find the coordinates of the images on the right (only need to change
    %x-coordinates);
    right_coords = positions;
    right_coords(:,1) = jth_right_image_Xcoords(right_coords(:,1), left_plane, ...
        right_plane, iterNum);    
    
    N = length(charges);
    Emat = zeros(N);
    
    %Define some constants:
    coulomb_constant = 8987551787.3681764; %in N m^2/C^2
    elementary_charge = 1.602176634E-19; %in C
    joules_to_eV = 6.242E+18;              
    
    %Accumulate interaction energy by considering every pair-wise
    %interaction between original partial charge and its two images in this
    %iteration
    for i = 1:N
        for j = 1:N
            
            %Find distances (still in nm)
            distL = sqrt(sum((positions(i,:) - left_coords(j,:)).^2));
            distR = sqrt(sum((positions(i,:) - right_coords(j,:)).^2));
            
            %Add both interactions energies, taking account of how the
            %image charge signs flip in each iteration (units are nm and e
            %still); also, crucial to divide by 2 to account for electric
            %field inside the electrodes being zero!!!
            Emat(i,j) = Emat(i,j) + (charges(i)*charges(j)*(-1)^iterNum)/(2*distL) + ...
                (charges(i)*charges(j)*(-1)^iterNum)/(2*distR);
            
        end
    end
    
    %Add in the Coulomb constant, convert charge from units of e to units
    %of Coulomb, and convert distance from nm to meters
    Emat = Emat * coulomb_constant * elementary_charge^2 * 1E9;
    
    %Convert from Joules to eV
    Emat = Emat * joules_to_eV;

end