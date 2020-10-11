function [smaller_HasMatch, larger_MatchID, reOrderSmaller, ...
    reOrderLarger, coords_SmallerSet] = match_up_atoms(...
    coords_SmallerSet, atomic_nums_SmallerSet, coords_LargerSet, ...
    atomic_nums_LargerSet,ToPlot)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Given two sets of similar atomic coordinates, 
    %we first try flipping the molecules with respect to each other
    %to determine the orientation in which they are the most similar. Then,
    %based on BOTH atom types AND atomic coordinates, we "match up" the
    %"equivalent" atoms from the two molecules and determine which atoms
    %are unique to each molecule. Returns information to indicate which
    %atoms match with each other. 
    %
    %~~~INPUTS~~~:
    %
    %coords_SmallerSet: 3-column array holding cartesian coordinates for
    %   each atom in the molecule with fewer atoms
    %
    %atomic_nums_SmallerSet: vector indicating the atomic number of each
    %   atom in the molecule with fewer atoms
    %
    %coords_LargerSet: 3-column array holding cartesian coordinates for
    %   each atom in the molecule with more atoms
    %
    %atomic_nums_LargerSet: vector indicating the atomic number of each
    %   atom in the molecule with more atoms
    %
    %ToPlot: logical variable; if true, a plot will be made showing which
    %   atoms were paired up with each other and which were determined to
    %   be "unique"
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    % 
    %smaller_HasMatch: a logical vector the same length as the number of
    %   atoms in the molecule with fewer atoms. Each position is true if
    %   and only if that atom was determined to have an "equivalent" or
    %   "matching" atom in the other molecule.
    %
    %larger_MatchID: a vector corresponding to each atom in the molecule
    %   with more atoms. If an atom was "matched" with an atom in the other
    %   molecule, then the position for that atom contains the row # of the
    %   atom in the smaller molecule that it was matched with. If not, then
    %   the position for that atom will be zero. 
    %
    %reOrderSmaller/reOrderLarger: vectors that can be used to re-order the
    %   coordinates of the smaller and larger molecules, respectively, to
    %   put the "matching" atoms first in the same order and the
    %   non-matching atoms at the end of each array.
    %
    %coords_SmallerSet: same as the input for this variable, except that
    %   if any flipping was done to make the two molecules align better
    %   that flipping will be reflected in these output coordinates.
    
    
    if nargin < 5
        ToPlot = false;
    end

    %First, figure out which way to orient the two molecules so that they
    %match up as good as possible; try 8 possible 180 degree flips along
    %x-, y-, and z- axes.
    %First we will make all 8 combinations of coordinates:
    smallerSet_testCoords = cell(8,1);
    for i = 1:8
        smallerSet_testCoords{i} = coords_SmallerSet;
    end
    for i = 5:8
        smallerSet_testCoords{i}(:,1) = -1*smallerSet_testCoords{i}(:,1);
    end
    for i = [3,4,7,8]
        smallerSet_testCoords{i}(:,2) = -1*smallerSet_testCoords{i}(:,2);
    end
    for i = [2,4,6,8]
        smallerSet_testCoords{i}(:,3) = -1*smallerSet_testCoords{i}(:,3);
    end    
    
    %Now we find the minimum distance between each set of test coordinates
    %and the other molecules
    minDists = zeros(8,1);
    for i = 1:8
        dists = pdist2(smallerSet_testCoords{i},coords_LargerSet);
        md = 0;
        for j = 1:size(dists,1)
            md = md + min(dists(j,:));
        end
        minDists(i) = md;
    end
    
    %Finally we select the test coordinates that were closest to the other
    %molecule
    [~,minIndex] = min(minDists);
    coords_SmallerSet = smallerSet_testCoords{minIndex};
        

    nSmall = size(coords_SmallerSet,1);
    nLarge = size(coords_LargerSet,1);
    
    smaller_HasMatch = false(nSmall,1);
    larger_MatchID = zeros(nLarge,1);
    
    dists = pdist2(coords_SmallerSet,coords_LargerSet);

    %Find smallest distance and the indicies corresponding to it
    [minDist,minRow] = min(dists);
    [minDist,minCol] = min(minDist);
    minRow = minRow(minCol);
    
    %Set threshhold to two times the median of the smallest distance
    %for each atom fromt the smaller set to the other molecule
    threshhold = 4*median(min(dists,[],2));
    
    while minDist < Inf
        
        if atomic_nums_SmallerSet(minRow) == atomic_nums_LargerSet(minCol) && ...
            minDist < threshhold
            smaller_HasMatch(minRow) = true;
            larger_MatchID(minCol) = minRow;
            dists(minRow,:) = Inf;
            dists(:,minCol) = Inf;
        else
            dists(minRow,minCol) = Inf;
        end
        
        %Re-calculate smallest distance and its indices
        [minDist,minRow] = min(dists);
        [minDist,minCol] = min(minDist);
        minRow = minRow(minCol);        
        
    end
    
    %Make vectors that show how to re-oder each list so that the conserved
    %atoms are put first, and all other atoms put last
    reOrderSmaller = zeros(nSmall,1);
    reOrderLarger = zeros(nLarge,1);
    
    %Make vector for re-ordering smaller list
    count_forward = 0;
    count_backward = nSmall+1;
    for i = 1:nSmall
        if smaller_HasMatch(i)
            count_forward = count_forward + 1;
            reOrderSmaller(count_forward) = i;
        else
            count_backward = count_backward - 1;
            reOrderSmaller(count_backward) = i;
        end
    end
    
    %Make vector for re-ordering larger list
    count_backward = nLarge+1;
    for i = 1:nLarge
        if larger_MatchID(i) > 0
            reOrderLarger(reOrderSmaller == larger_MatchID(i)) = i;
        else
            count_backward = count_backward - 1;
            reOrderLarger(count_backward) = i;
        end
    end

    if ToPlot
        figure();
        hold on;
        for i = 1:nSmall
            if smaller_HasMatch(i)
                text(coords_SmallerSet(i,1),coords_SmallerSet(i,2),...
                    coords_SmallerSet(i,3),num2str(i),'Color','b');
            else
                text(coords_SmallerSet(i,1),coords_SmallerSet(i,2),...
                    coords_SmallerSet(i,3),'0','Color','b');                
            end
        end
        for i = 1:nLarge
            text(coords_LargerSet(i,1),coords_LargerSet(i,2),...
                coords_LargerSet(i,3),num2str(larger_MatchID(i)),'Color','r');
        end
        
        minL = min(min(coords_LargerSet));
        maxL = max(max(coords_LargerSet));
        xlim([minL maxL]);
        ylim([minL maxL]);
        zlim([minL maxL]);
        
        xlabel('X (nm)');
        ylabel('Y (nm)');
    end
    
end