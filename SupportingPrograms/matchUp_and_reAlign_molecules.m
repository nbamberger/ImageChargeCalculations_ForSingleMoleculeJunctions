function [coords1, coords2, an1, an2, ch1, ch2, r1, r2, nConserved] = ...
    matchUp_and_reAlign_molecules(coords1, coords2, an1, an2, ch1, ch2, ...
    ToPlot,DoNotMatch_IDList)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Similar to match_up_atoms, but instead of just
    %figuring out which atoms correspond to which other atoms, this 
    %function also re-orders and potentially flips the coordinates (and
    %atomic numbers and charges) so that the two lists of atoms match up as
    %well as possible, with any non-matching atoms at the end. For this
    %function, it doesn't matter whether molecule #1 or molecule #2 has
    %more atoms in it. 
    %
    %~~~INPUTS~~~:
    %
    %coords1/coords2: 3-column array containing the cartesian coordinates
    %   for each atom in molecule #1/#2
    %
    %an1/an2: vector listing the atomic number for each atom in molecule
    %   #1/#2
    %
    %ch1/ch2: vector listing the partial charge assigned to each atom in
    %   molecule #1/#2
    %
    %ToPlot: logical variable; if true, a plot will be made showing which
    %   atoms were paired up with each other and which were determined to
    %   be "unique"
    %
    %DoNotMatch_IDList: an optional vector containing the ID #s of atoms
    %   (from the smaller coordinate set) that should NOT be matched with
    %   the larger molecule, no matter if they are close enough or not
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %coords1/coords2: same as the input coordinate arrays, but with the
    %   rows (atoms) re-ordered to put the "matching" atoms between the two
    %   molecules at the beginning and in the same order, and the
    %   non-matching atoms at the end of each array
    %
    %an1/an2: again, same as input but with atoms re-ordered to put
    %   matching atoms at the beginning in the same order
    %
    %ch1/ch2: same as input but with atoms re-ordered, etc.
    %
    %r1/r2: the vectors that were used to re-order the inputs to produce
    %   the outputs; so, for example, in the output an1 = an1(r1) and an2 =
    %   an2(r2)
    %
    %nConserved: the number of atoms that are "conserved" between the two
    %   molecules, i.e. how many atoms were matched up as opposed to being
    %   considered unique to one or the other molecule.
    
    
    if nargin < 7
        ToPlot = false;
    end
    if nargin < 8
        DoNotMatch_IDList = [];
    end

    %Run match-up atoms, with the order of inputs depending on which list
    %of atoms is smaller.  For the smaller set, replace the coordinates
    %with the output fro match-up atoms so that any flips are accounted for
    if length(ch1) > length(ch2)
        [conserved,~,r2,r1,coords2] = match_up_atoms(coords2,an2,coords1,...
            an1,ToPlot,DoNotMatch_IDList);
    else
        [conserved,~,r1,r2,coords1] = match_up_atoms(coords1,an1,coords2,...
            an2,ToPlot,DoNotMatch_IDList);
    end
    
    %Re-order all inputs to get final outputs:
    coords1 = coords1(r1,:);
    coords2 = coords2(r2,:);
    an1 = an1(r1);
    an2 = an2(r2);
    ch1 = ch1(r1);
    ch2 = ch2(r2);
    
    %Calculate # of conserved atoms between the two molecules
    nConserved = sum(conserved);

end