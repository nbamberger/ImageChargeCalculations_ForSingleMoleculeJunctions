function [nConserved, Emat1, Emat2, ch1, ch2, pos1, pos2, an1, an2] = ...
    prepare_2mol_imagecharge_comparison(charges1, charges2, positions1, ...
    positions2, atomic_nums1, atomic_nums2, gap_size, minEnChange,...
    DoNotMatch_IDList)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Given coordinate, atom type, and partial charge
    %info for two similar molecules, this function carries out an image
    %charge re-normalization calculation for both molecules. It also
    %matches up the two atoms to determine which atoms are conserved
    %between the two and which are unique to one molecule or the other.
    %Finally, it re-orders the renormalization outputs to list the
    %conserved atoms first. This is useful for later comparisons of the two
    %molecules' renormalizations. 
    %
    %~~~INPUTS~~~:
    %
    %charges1/charges2: vector listing the partial charge assigned to each
    %   atom in molecule #1/#2 in units of e
    %
    %positions1/positions2: 3-column array listing cartesian coordinates
    %   for each atom in molecule #1/#2 in units of nm
    %
    %atomic_nums2/atomic_nums2: vector listing the atomic number for each
    %   atom in molecule #1/#2
    %
    %gap_size: the size of the gap between image planes, in units of nm
    %
    %minEnChange: higher order charge images will continue to be added
    %   until the total cumulative interaction energy varies by less than
    %   this value (in units of eV)
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %nConserved: the number of atoms conserved between the two molecules
    %
    %Emat1/Emat2: atomic interaction matrices for molecule 1/2, in which
    %   the (i,j)th element contains the image charge interaction energy
    %   between atom i and all images of atom j. The atoms have been
    %   re-ordered to put the atoms conserved between the two molecules
    %   first.
    %
    %ch1/ch2: same as input charges but re-ordered to put conserved atoms
    %   first
    %
    %pos1/pos2: same as coordinate inputs but re-ordered to put conserved
    %   atoms first
    %
    %an1/an2: same as atomic number inputs but re-ordered to put conserved
    %   atoms first
    %
    %DoNotMatch_IDList: an optional vector containing the ID #s of atoms
    %   (from the smaller coordinate set) that should NOT be matched with
    %   the larger molecule, no matter if they are close enough or not
    
    
    if nargin < 9
        DoNotMatch_IDList = [];
    end
    if nargin < 8
        minEnChange = 0.002;
    end
    if nargin < 7
        gap_size = 2.05*max(abs(positions2(:,1)));
    end

    %Calculate contributions to total renormalization energy from each
    %pairwise atomic combination, for each molecule
    [~,Emat1] = Detailed_TwoElectrode_Renormalization(charges1,positions1,...
        gap_size,atomic_nums1,'AllSeparate',false,minEnChange,false);
    [~,Emat2] = Detailed_TwoElectrode_Renormalization(charges2,positions2,...
        gap_size,atomic_nums2,'AllSeparate',false,minEnChange,false);    
    
    %Match up atoms that are "the same" between the two molecules and
    %re-order the atoms so that the conserved atoms all come first
    [pos1,pos2,an1,an2,ch1,ch2,r1,r2,nConserved] = matchUp_and_reAlign_molecules(...
        positions1,positions2,atomic_nums1,atomic_nums2,charges1,charges2,...
        false,DoNotMatch_IDList);
        
    %Re-order matrices to put conserved coordinates first
    Emat1 = Emat1(r1,:);
    Emat1 = Emat1(:,r1);
    Emat2 = Emat2(r2,:);
    Emat2 = Emat2(:,r2);

end