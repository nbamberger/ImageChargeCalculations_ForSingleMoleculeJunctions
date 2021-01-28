function DisplaySystem_SphericalElectrodes(charges, positions, ...
    left_center, right_center, left_radius, right_radius, nIter)
    %Copyright 2021 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Makes a figure to display the charge system in
    %the case of spherical electrodes, showing the electrodes, the original
    %charges, and image charges inside those electrodes if requested.
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %positions: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %left_center/right_center: the x-coordinate center (in nm) of the 
    %   left/right spherical electrode; the y- and z-coordinates of both
    %   centers are assumed to be zero, i.e. both centers lie on the x-axis
    %
    %left_radius/right_radius: the radius (in nm) of the left/right
    %   spherical electrode
    %
    %nIter: how many iterations of image charges to display; can be set to
    %   zero   
    
    
    size_scalar = 50;
    
    figure();
    hold on;

    %Colors for positive and negative charges
    cols = [0 0 1; 0 0 0; 1 0 0];
    
    %Plot the original charges
    scatter(positions(:,1),positions(:,2),abs(charges)*size_scalar,...
        cols(sign(charges)+2,:),'filled');
    xlabel('X (nm)');
    ylabel('Y (nm)');
    
    %Display the spherical electrodes
    gold = [0.7 0.7 0.25];
    t = linspace(0,2*pi(),100);
    fill([cos(t)*left_radius + left_center,left_radius+left_center],...
        [sin(t)*left_radius,0],gold);
    fill([cos(t)*right_radius + right_center,right_radius+right_center],...
        [sin(t)*right_radius,0],gold);    
    
    %Show image charges, if requested
    if nIter > 0
        
        [left_ch, left_co] = getSphere_ImageCharges(charges, positions, ...
            left_center, left_radius);
        [right_ch, right_co] = getSphere_ImageCharges(charges, positions, ...
            right_center, right_radius);   
        
        scatter(left_co(:,1),left_co(:,2),abs(right_ch)*size_scalar,...
            cols(sign(left_ch)+2,:),'filled');
        scatter(right_co(:,1),right_co(:,2),abs(right_ch)*size_scalar,...
            cols(sign(right_ch)+2,:),'filled');
        
        for i = 2:nIter
            
            %Get the next order of image charges
            left_ch_old = left_ch;
            left_co_old = left_co;      
            [left_ch, left_co] = getSphere_ImageCharges(right_ch, right_co, ...
                left_center, left_radius);
            [right_ch, right_co] = getSphere_ImageCharges(left_ch_old, left_co_old, ...
                right_center, right_radius); 
            
            scatter(left_co(:,1),left_co(:,2),abs(right_ch)*size_scalar,...
                cols(sign(left_ch)+2,:),'filled');
            scatter(right_co(:,1),right_co(:,2),abs(right_ch)*size_scalar,...
                cols(sign(right_ch)+2,:),'filled');            
            
        end
        
    end
    
    %Make sure the axes are plotted on equivalent scales:
    a = gca;
    xL = a.XLim;
    yL = a.YLim;
    mm = max(abs([xL,yL]));
    xlim([-mm mm]);
    ylim([-mm mm]);
    a.Units = 'inches';
    minlength = min(a.Position(3:4));
    a.Position(3) = minlength;
    a.Position(4) = minlength;

end