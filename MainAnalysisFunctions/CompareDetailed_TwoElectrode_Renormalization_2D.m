function [E1,E2,Ediff] = ...
    CompareDetailed_TwoElectrode_Renormalization_2D(charges1, charges2, ...
    coords1, coords2, atomic_nums1, atomic_nums2,gap_size, minEnChange,...
    CombinedSwitch, ConservedBuckets, CombineAcrossDiag)
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
    %and divvies up this energy between each unique PAIR of atoms (e.g.,
    %the renormalization energy due to hydrogen1 interacting with itself,
    %due to hydrogen1 interacting with hydrogen2, etc.). Finally, a 2D
    %matrix is plotted to show how the renormalization energy of different
    %types of interactions changes going from molecule #1 to #2. The
    %interactions of the conserved atoms can be grouped in different ways.
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
    %
    %CombinedSwitch: logical variable; if true, the effects or removing any
    %   unique atoms from molecule 1 and of adding any unique atoms to
    %   moleucle 2 will be combined into a single bucket, if false the two
    %   will be listed separately.
    %
    %ConservedBuckets: the name of the way the user wishes to group the
    %   conserved atoms for comparison purposes. Options are "AllTogether"
    %   (single bucket for all conserved atoms), "HydrogensSeparate" (one
    %   bucket for conserved hydrogens, one bucket for other conserved
    %   atoms), "ByAtomType" (one bucket for each conserved atomic
    %   species), and "AllSeparate" (one bucket for each conserved atom).
    %
    %CombineAcrossDiag: logical variable; whether or not to combine buckets
    %   across the diagonal (i.e. add together the interactions of
    %   [hydrogen1 with images of hydrogen2] and [hydrogen2 with images of
    %   hydrogen1] and only list the sum once). Due to the inherent
    %   symmetry of the image charge forumla, the interaction matrix will
    %   ALWAYS be symmetric across the diagonal, and so this combination
    %   does not lose any information.
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %E1/E2: the matrix of renormalization energies for molecule 1/2, 
    %   grouped into buckets as specified by the user
    %
    %Ediff: the matrix representing the changes in renormalization
    %   energies; simply E2 - E1
    
    
    %Default inputs
    if nargin < 11
        CombineAcrossDiag = true;
    end
    if nargin < 10
        ConservedBuckets = 'AllTogether';
    end
    if nargin < 9
        CombinedSwitch = true;
    end
    if nargin < 8
        minEnChange = 0.002;
    end
    if nargin < 7
        gap_size = 2;
    end
    
    %Calcualte the total pairwise matrices for each molecule and re-order
    %them to put all conserved atoms first
    [nConserved, Emat1, Emat2, ~, ~, ~, ~, an1, an2] = ...
        prepare_2mol_imagecharge_comparison(charges1, charges2, coords1, ...
        coords2, atomic_nums1, atomic_nums2, gap_size, minEnChange);
    n1 = length(an1);
    n2 = length(an2);
    
    if CombinedSwitch
        nBucketsSub = 1;
    else
        nBucketsSub = 2;
    end
    
    %Load in atomic symbol dictionary
    atomic_symbol_dictionary = importdata('atomic_symbol_dictionary.mat');
    
    %Now we need to assign each conserved atom to a bucket #
    if strcmp(ConservedBuckets,'AllTogether')
        
        nBuckets = 1 + nBucketsSub;
        bucket_IDs = ones(nConserved,1);
        
        axis_labels = cell(nBuckets,1);
        axis_labels{1} = 'Conserved Atoms';
        
    elseif strcmp(ConservedBuckets,'ByAtomType')
        
        [atom_types,~,bucket_IDs] = unique(an1(1:nConserved));
        
        nAtoms = length(atom_types);
        nBuckets = nAtoms + nBucketsSub;
        
        axis_labels = cell(nBuckets,1);
        axis_labels(1:nAtoms) = atomic_symbol_dictionary(atom_types);
        
    elseif strcmp(ConservedBuckets,'HydrogensSeparate')
        
        nBuckets = 2 + nBucketsSub;
        
        bucket_IDs = zeros(nConserved,1);
        for i = 1:nConserved
            if an1(i) == 1
                bucket_IDs(i) = 1;
            else
                bucket_IDs(i) = 2;
            end
        end
        
        axis_labels = cell(nBuckets,1);
        axis_labels{1} = 'Conserved Hs';
        axis_labels{2} = 'Other Conserved Atoms';
        
    elseif strcmp(ConservedBuckets,'AllSeparate')
        
        nBuckets = nConserved + nBucketsSub;
        bucket_IDs = (1:nConserved)';
        
        axis_labels = cell(nBuckets,1);
        axis_labels(1:nConserved) = atomic_symbol_dictionary(an1(1:nConserved));
        
    else
        error('Unrecognized input for ConservedBuckets; valid inputs are "AllTogether", "ByAtomType", "HydrogensSeparate", "AllSeparate"');
    end
    
    %Now we can assign ALL atoms to buckets
    bucket_IDs_1 = zeros(n1,1);
    bucket_IDs_2 = zeros(n2,1);
    bucket_IDs_1(1:nConserved) = bucket_IDs;
    bucket_IDs_2(1:nConserved) = bucket_IDs;
    %Substituent 1 atoms will go in the second to last bucket if the
    %substituents are not being combined, and the last bucket if they are
    if CombinedSwitch
        bucket_IDs_1(nConserved+1:n1) = nBuckets;
        axis_labels{nBuckets} = 'Substituent Switch';
    else
        bucket_IDs_1(nConserved+1:n1) = nBuckets-1;
        axis_labels{nBuckets} = 'Mol. 2 Substituent';
        axis_labels{nBuckets - 1} = 'Mol. 1 Substituent';
    end
    %Substituent 2 atoms will go in the last bucket no matter what
    bucket_IDs_2(nConserved+1:n2) = nBuckets;
    
    
    %Now we can create the summary matrices by grouping values into their
    %appropriate buckets
    E1 = zeros(nBuckets);
    E2 = zeros(nBuckets);
    
    %Put all atoms into correct buckets
    for i = 1:n1
        for j = 1:n1
            E1(bucket_IDs_1(i),bucket_IDs_1(j)) = ...
                E1(bucket_IDs_1(i),bucket_IDs_1(j)) + Emat1(i,j);          
        end
    end
    for i = 1:n2
        for j = 1:n2
            E2(bucket_IDs_2(i),bucket_IDs_2(j)) = ...
                E2(bucket_IDs_2(i),bucket_IDs_2(j)) + Emat2(i,j);            
        end
    end
    
    %Calculate difference matrix
    Ediff = E2 - E1;
    
    %Get orange-teal colormap
    cmap = importdata('orange_cyan_cmap.mat');
    
    %Make plot for molecule 2 minus moleclue 1
    hotcold_heatmap_of_matrix(Ediff,CombineAcrossDiag,cmap);
    xticklabels(axis_labels);
    yticklabels(axis_labels);
    title('Molecule 2 - Molecule1');
    xlabel('Change interaction of partial charges on...');
    ylabel('...with image charges of partial charges from:')   
    
    %Label color bar
    f = gcf;
    f.Children(1).Label.String = 'Interaction Energy (eV)';
    f.Children(1).Label.Rotation = 270;
    f.Children(1).Label.VerticalAlignment = 'bottom';

end