function DifferenceInPartialCharge_PerAtom(charges1, charges2, ...
    coords1, coords2, atomic_nums1,atomic_nums2,DoNotMatch_IDList)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Given two very similar molecules (e.g. only
    %differing at a single substituent), this function first compares the
    %two molecules to determine which atoms are "conserved"--i.e. belong to
    %parts of the molecule that are the same as each other--and which atoms
    %are unique to either molecule #1 or molecule #2. This function then
    %makes a plot showing the change in partial charge on each conserved
    %atom going from molecule 1 to 2. The non-conserved atoms are not
    %plotted, but their average position is indicated with an x. For the
    %circles, the area is proportional the magnitude of change in charge,
    %and orange indicates a positive change and teal a negative change.
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
    %DoNotMatch_IDList: an optional vector containing the ID #s of atoms
    %   (from the smaller coordinate set) that should NOT be matched with
    %   the larger molecule, no matter if they are close enough or not
    
    
    if nargin < 7
        DoNotMatch_IDList = [];
    end
    
    %Match up atoms in order to put conserved atoms first
    [nConserved, ~, ~, ch1, ch2, pos1, pos2, ~, ~] = ...
        prepare_2mol_imagecharge_comparison(charges1, charges2, coords1, ...
        coords2, atomic_nums1, atomic_nums2, 10, 100,DoNotMatch_IDList);
    n2 = length(ch2);
    n1 = length(ch1);

    %Add a plot to show the change in partial charge for each conserved
    %atom; make sure to use square root appropriately so that circle AREA
    %is proportional to change in partial charge
    figure();
    hold on;
    circle_scalar = 50;
    cols = [0 1 1;0 0 0; 1 0.5 0];
    for i = 1:nConserved
        ch_diff = ch2(i) - ch1(i);
        col_index = sign(ch_diff) + 2;
        plot(pos1(i,1),pos1(i,2),'o','Color',cols(col_index,:),...
            'MarkerFaceColor',cols(col_index,:),'MarkerSize',abs(sqrt(ch_diff))*circle_scalar);
    end
    title('Change in Partial Charges on Conserved Atoms');
    xlabel('X (nm)');
    ylabel('Y (nm)');
    
    %Add an OPEN circle for the substituent (at the average position of the
    %substituent in molecule 2)
    pos = mean(pos2(nConserved+1:n2,:),1);
    ch_diff = sum(ch2(nConserved+1:n2)) - sum(ch1(nConserved+1:n1));
    col_index = sign(ch_diff) + 2;
    plot(pos(1),pos(2),'o','Color',cols(col_index,:),'MarkerSize',...
        abs(sqrt(ch_diff))*circle_scalar);
    
    %Adjust bounds
    mm = max(max(abs(pos1(:,1:2))))*1.1;
    xlim([-mm mm]);
    ylim([-mm mm]);

end