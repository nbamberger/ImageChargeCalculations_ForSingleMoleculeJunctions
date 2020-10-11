function DifferenceInRenormEnergy_PerAtom(charges1, charges2, ...
    coords1, coords2, atomic_nums1,atomic_nums2,gap_size,minEnChange)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Given two very similar molecules (e.g. only
    %differing at a single substituent), this function first compares the
    %two molecules to determine which atoms are "conserved"--i.e. belong to
    %parts of the molecule that are the same as each other--and which atoms
    %are unique to either molecule #1 or molecule #2. Next, this function
    %calculates the image charge renormalization energy for each molecule
    %and divvies up this energy between all the atoms. This function then
    %makes a plot showing the change in renormalization energy on each
    %conserved atom going from molecule 1 to 2. The non-conserved atoms are
    %combined into a single open circle plotted at their average position.
    %For the circles, the area is proportional the magnitude of change in
    %renormalization energy, and red/blue indicates a positive/negative
    %change.
    %
    %~~~INPUTS~~~:
    %
    %charges1/charges2: a vector listing the partial charge (in units of e)
    %   assigned to each atom in molecule 1/2
    %
    %coords1/coords2: a 3-column array listing the x-, y-, and
    %   z-coordinates (in units of nm) of each atom in molecule 1/2
    %
    %atomic_nums1/atomic_nums2: a vector listing the atomic number of each
    %   atom in molecule 1/2. This information is needed to propertly match
    %   up which atoms are "equivalent" or not between the two molecules
    %
    %gap_size: the distance between the two electrodes (or, to be exact,
    %   the image planes of the two electrodes) in units of nm
    %
    %minEnChange: a minimum energy change in units of eV. Higher order
    %   image charges will be added until the change in energy is less than
    %   this value
    

    %Default inputs
    if nargin < 8
        minEnChange = 0.002;
    end
    if nargin < 7
        gap_size = 2;
    end
    
    %Calcualte the total pairwise matrices for each molecule and re-order
    %them to put all conserved atoms first
    [nConserved, Emat1, Emat2, ~, ~, pos1, pos2, an1, an2] = ...
        prepare_2mol_imagecharge_comparison(charges1, charges2, coords1, ...
        coords2, atomic_nums1, atomic_nums2, gap_size, minEnChange);
    n1 = length(an1);
    n2 = length(an2);
    
    %Marginalize over matrices to get renormalization energy "assigned" to
    %each atom
    Esum1 = sum(Emat1);
    Esum2 = sum(Emat2);

    %Add a plot to show the change in renormalization energy for each
    %conserved atom
    figure();
    hold on;
    circle_scalar = 100;
    cols = [0 1 1;0 0 0; 1 0.5 0];
    for i = 1:nConserved
        en_diff = Esum2(i) - Esum1(i);
        col_index = sign(en_diff) + 2;
        plot(pos1(i,1),pos1(i,2),'o','Color',cols(col_index,:),...
            'MarkerFaceColor',cols(col_index,:),'MarkerSize',sqrt(abs(en_diff))*circle_scalar);
    end
    
    %Add an open circle to show the change in energy due to the substituent
    %switch
    en_diff = sum(Esum2(nConserved+1:n2)) - sum(Esum1(nConserved+1:n1));
    col_index = sign(en_diff) + 2;
    pos = mean(pos2(nConserved+1:n2,:),1); %Average position of atoms in second substituent   
    plot(pos(1),pos(2),'o','Color',cols(col_index,:),'MarkerSize',...
        max(sqrt(abs(en_diff))*circle_scalar,1));
    
    title('Change in Renormalization Energy for Each Atom');
    xlabel('X (nm)');
    ylabel('Y (nm)');
    mm = max(max(abs(pos1(:,1:2))))*1.1;
    xlim([-mm mm]);
    ylim([-mm mm]);

end