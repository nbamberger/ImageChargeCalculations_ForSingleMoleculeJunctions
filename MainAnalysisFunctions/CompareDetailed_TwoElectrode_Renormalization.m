function [Ecomparison,labels,conserved] = ...
    CompareDetailed_TwoElectrode_Renormalization(charges1, charges2, ...
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
    %and divvies up this energy between all the atoms. Finally, plots are
    %made to show how the renormalization energy of the conserved atoms
    %changes going from molecule #1 to #2. 
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
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %Ecomparison: a two-column array listing the image charge
    %   renormalization energy assigned to each atom in molecule 1/2 in the
    %   1st/2nd column, respectively. Conserved atoms are listed first
    %   (and, since they are conserved, values for both molecules are in
    %   the same row), followed by the atoms unique to molecule 1 and
    %   finally the atoms unique to molecule 2.
    %
    %labels: a one-column celly array listing the atomic symbol
    %   corresponding to each row of Ecomparison
    %
    %conserved: a logical vector listing whether each row of Ecomparison is
    %   for a conserved atom (true) or an atom unique to one of the
    %   molecules (false)
    
    
    %Default inputs
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

    %Marginalize over matrices to get renormalization energy "assigned" to
    %each atom    
    Esum1 = sum(Emat1);
    Esum2 = sum(Emat2);

    %The number of different points to plot is the sum of the number of
    %conserved atoms in both molecules plus the atoms only in molecule 1 or
    %only in molecule 2
    nPlot = nConserved + (n1 - nConserved) + (n2 - nConserved);
    
    %Get summary comparison data to return to user
    Ecomparison = zeros(nPlot,2);
    Ecomparison(1:n1,1) = Esum1;
    Ecomparison(1:nConserved,2) = Esum2(1:nConserved);
    Ecomparison(n1+1:nPlot,2) = Esum2(nConserved+1:n2);
    
    conserved = false(nPlot,1);
    conserved(1:nConserved) = true;
    
    %Make plot and add conserved atoms
    figure();
    hold on;
    plot((1:nConserved),Esum1(1:nConserved),'-o','Color','b');
    plot((1:nConserved),Esum2(1:nConserved),'-o','Color','r');
    xlabel('Atom');
    ylabel('Renormalization Energy (eV)');
    legend({'Molecule 1','Molecule 2'},'Location','best','AutoUpdate','off');

    %Add non-conserved atoms
    counter = nConserved;
    for i = nConserved+1:n1
        counter = counter + 1;
        plot(counter, Esum1(i),'o','Color','b');
    end
    for i = nConserved+1:n2
        counter = counter + 1;
        plot(counter, Esum2(i),'o','Color','r');
    end     
    
    %Generate atom labels
    labels = zeros(nPlot,1);
    labels(1:n1) = an1;
    labels(n1+1:nPlot) = an2(nConserved+1:n2);   
    atomic_symbol_dictionary = importdata('atomic_symbol_dictionary.mat');
    labels = atomic_symbol_dictionary(labels);

    xticks((1:nPlot));
    xticklabels(labels);
    
    hold on;
    plot([1,nPlot],[0,0],'--','Color',[0.5 0.5 0.5]);
    a = gca;
    a.XGrid = 'on';
    
    %Add plot to show the difference in energy for each atom
    figure();
    plot((1:nPlot),Ecomparison(:,2) - Ecomparison(:,1),'-o','Color',[0 0 0]);
    xticks((1:nPlot));
    xticklabels(labels);
    xlabel('Atom');
    ylabel('\Delta in Renorm. Energy, Molecule 2 - Molecule 1 (eV)');
    
    hold on;
    plot([1,nPlot],[0,0],'--','Color',[0.5 0.5 0.5]);
    mm = max(abs(Ecomparison(:,2) - Ecomparison(:,1)));
    ylim([-mm mm]);
    plot([nConserved+0.5 nConserved+0.5],[-mm mm],':','Color','g','LineWidth',2);  

end