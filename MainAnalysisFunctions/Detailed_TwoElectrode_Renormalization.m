function [TotalEnergy, EnergyMatrix] = ... 
    Detailed_TwoElectrode_Renormalization(charges, positions, gap_size, ...
    atom_nums, BucketType, CombineAcrossDiag, minEnChange, ToPlot)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Similar to SingleGapSize_Image
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %positions: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %gap_size: the gap between the two electrodes, in units of nm.  Each
    %   electrode is assumed to be in the yz-plane and the two are assumed
    %   to be symmetric about zero.  Thus, for example, a gap size of 1
    %   means that the two electrodes are located at -0.5 and 0.5 nm.  
    %
    %atom_names: a vector listing the atomic number for each atom; can be
    %   left out
    %
    %BucketType: How the atoms should be grouped for the matrix; options
    %   are 'AllSeparate' or 'ByAtom'
    %
    %CombineAcrossDiag: Logical variable; whether or not to combine
    %   energies across the diagonal of the matrix
    %
    %minEnChange: This is the convergence limit, in eV, for the total image
    %   charge energy. In other words, the calculation will keep adding
    %   more charge images until the change in total energy is less than
    %   this value
    %
    %ToPlot: logical variable, whether or not to create visible plots
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %    
    %TotalEnergy: the total cumulative image charge renormalization energy
    %   after the final iteration
    %
    %EnergyMatrix: A matrix in which the (i,j)th element contains the
    %   interaction energy (after the final iteration) between the ith
    %   partial charge and all images of the jth partial charge.  
    
    
    %Default inputs
    if nargin < 8
        ToPlot = true;
    end
    if nargin < 7
        minEnChange = 0.002;
    end
    if nargin < 6
        CombineAcrossDiag = false;
    end
    if nargin < 5
        BucketType = 'AllSeparate';
    end
    if nargin < 4
        %If atom names are missing, cannot group by atom!
        atom_nums = [];
        BucketType = 'AllSeparate';
    end
    
    if ~any(strcmp(BucketType,{'AllSeparate','ByAtom'}))
        error('Unrecognized bucket type; options are "AllSeparate" or "ByAtom"');
    end
    
    if ~isempty(atom_nums)
        %Load in atomic symbol dictionary
        atomic_symbol_dictionary = importdata('atomic_symbol_dictionary.mat');
    end
    
    nAtoms = length(charges);
    EnergyMatrix = zeros(nAtoms);
        
    %X-coordinates of the two image planes
    image_plane1 = -1*gap_size/2;
    image_plane2 = gap_size/2;
    
    %Display a warning if any of the partial charges are on the wrong side
    %of the image planes
    if min(positions(:,1)) < image_plane1 || max(positions(:,1)) > ...
            image_plane2
        warning('Atom(s) are inside electrodes!!');
    end
    
    %Keep adding more images until TOTAL energy converged to within limit
    maxIter = 3000;
    enChange = Inf;
    i = 0;
    while abs(enChange) >= minEnChange && i < maxIter
        i = i + 1;
        matChange = interaction_energy_with_jth_images_MATRIX(charges, ...
            positions,image_plane1,image_plane2,i);
        enChange = sum(sum(matChange));
        EnergyMatrix = EnergyMatrix + matChange;
    end
    
    %Alert user if conergence was not reached
    if abs(enChange) > minEnChange
        warning('Reached maximum # of iterations without meeting convergence criteria');
        disp(strcat('Threshhold:', {' '}, num2str(minEnChange), ' Last Change:', ...
            {' '}, num2str(enChange)));
    end
    
    %Add one more half-change, because we know that the the changes at each
    %step oscillate up and down
    EnergyMatrix = EnergyMatrix + 0.5 * interaction_energy_with_jth_images_MATRIX(...
        charges,positions,image_plane1,image_plane2,i+1);

    TotalEnergy = sum(sum(EnergyMatrix));
    
    %Make new matrix grouped by atom, if that is requested
    if strcmp(BucketType,'ByAtom')
        %Get a list of unique atoms and figure out which atoms belong in
        %which bin
        [PlotNames, ~, unique_IDs] = unique(atom_nums);
        nUnique = length(PlotNames);
        
        %Create the new matrix for plotting by assigning each atom pair to
        %the correct bin
        PlotMatrix = zeros(nUnique);
        for i = 1:nAtoms
            for j = 1:nAtoms
                PlotMatrix(unique_IDs(i),unique_IDs(j)) = ...
                    PlotMatrix(unique_IDs(i),unique_IDs(j)) + EnergyMatrix(i,j);
            end
        end
    else
        PlotMatrix = EnergyMatrix;
        PlotNames = atom_nums;
    end
    
    if ToPlot
        %Re-arrange data for plotting
        ColumnData = zeros(nAtoms^2,3);
        counter = 0;
        for i = 1:nAtoms
            for j = 1:nAtoms
                counter = counter + 1;
                ColumnData(counter,1) = i;
                ColumnData(counter,2) = j;
                ColumnData(counter,3) = EnergyMatrix(i,j);
            end
        end

        %Make figure to show matrix
        hotcold_heatmap_of_matrix(PlotMatrix,CombineAcrossDiag);
        xlabel('...Original Partial Charge:');
        ylabel('...And Images of Partial Charge:');
        title('Renormalization energy from interation of...')
        
        %Label color bar
        f = gcf;
        f.Children(1).Label.String = 'Interaction Energy (eV)';
        f.Children(1).Label.Rotation = 270;
        f.Children(1).Label.VerticalAlignment = 'bottom';

        if ~isempty(atom_nums)
            xticklabels(atomic_symbol_dictionary(PlotNames));
            yticklabels(atomic_symbol_dictionary(PlotNames));
            ytickangle(90);
        end

        %Calculate atom contributions (I *think* that the EnergyMatrix will
        %always be symmetric, so it doesn't matter how we divide up)
        atom_contributions = sum(EnergyMatrix);

        %Make figure to summarize by original atom
        figure();
        hold on;
        plot((1:nAtoms),atom_contributions,'-o');
        if ~isempty(atom_nums)
            xticks((1:nAtoms));
            xticklabels(atomic_symbol_dictionary(atom_nums));
        end    
        xlabel('Original Atom');
        ylabel('Contribution to Renormalization (eV)');
        plot([0,nAtoms],[0 0],'--','Color',[0.5 0.5 0.5]);
        a = gca;
        a.XGrid = 'on';

        %Make a plot showing the geometry of the system
        Display_ChargeElectrode_System(charges, positions, -(gap_size)/2, ...
            (gap_size)/2,0);
    end

end