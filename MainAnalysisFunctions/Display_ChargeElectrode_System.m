function Display_ChargeElectrode_System(charges, positions, left_plane, ...
    right_plane, nImageIterations)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Creates a plot showing where the partial 
    %charges, image charges, and image planes are located
    %
    %~~~INPUTS~~~:
    %
    %charges: vector of partial charges assigned to each atom, in units of
    %   e
    %
    %positions: 3-column matrix of x-, y-, and z- coordinates for each
    %   atom, in units of nm
    %
    %left_plane/right_plane: the y-coordinate of the left/right planar
    %   electrode, respectively (assumed to be in the yz-plane), in units
    %   of nm
    %
    %nImageIterations: the number of image charge iterations to display
    
    
    figure();
    hold on;
    N = length(charges);
    size_scalar = 75;
    
    %Determine colors for each charge (blue for negative, red for positive)
    C = [0 0 1; 1 0 0; 0 0 0];
    Cinv = [1 0 0; 0 0 1; 0 0 0];
    ch_signs = ones(N,1);
    ch_signs(charges > 0) = 2;
    ch_signs(charges == 0) = 3;
    
    %Make circle sizes proportional to absolute charge, but also set a
    %minimum size; use square root so that the AREA of the circle is the
    %thing proportional to the charge
    circ_sizes = sqrt(abs(charges))*size_scalar;
    circ_sizes = max(circ_sizes,1);
    
    %Plot the original charges:
    scatter(positions(:,1),positions(:,2),circ_sizes,C(ch_signs,:),'filled');
    boundsX = [min(positions(:,1)) max(positions(:,1))];
    boundsY = [min(positions(:,2)) max(positions(:,2))];
    
    %Plot the image charges:
    inversion = true;
    for i = 1:nImageIterations
        %Calculate positions of images on left and right
        coordsL = positions;
        coordsL(:,1) = jth_left_image_Xcoords(coordsL(:,1),left_plane,...
            right_plane, i);
        
        coordsR = positions;
        coordsR(:,1) = jth_right_image_Xcoords(coordsR(:,1),left_plane,...
            right_plane, i);        
        
        %Plot the images using appropriate color scheme
        if inversion
            scatter(coordsL(:,1),coordsL(:,2),circ_sizes,Cinv(ch_signs,:),'filled');
            scatter(coordsR(:,1),coordsR(:,2),circ_sizes,Cinv(ch_signs,:),'filled');
        else
            scatter(coordsL(:,1),coordsL(:,2),circ_sizes,C(ch_signs,:),'filled');
            scatter(coordsR(:,1),coordsR(:,2),circ_sizes,C(ch_signs,:),'filled');            
        end
        
        %Switch charge inversion for next iteration
        inversion = ~inversion;
        
        boundsX = [min(coordsL(:,1)) max(coordsR(:,1))];
        boundsY = [min(coordsL(:,2)) max(coordsR(:,2))];
    end

    %Adjust bounds
    newMax = 1.1*max([abs(boundsX), abs(boundsY), left_plane, right_plane]);
    bounds = [-newMax newMax];   
    xlim(bounds);
    ylim(bounds);
    
    %Plot the image planes
    plot([right_plane, right_plane], bounds, ...
        'Color', [0 0 0],'LineWidth',1.5);
    plot([left_plane, left_plane], bounds, ...
        'Color', [0 0 0],'LineWidth',1.5);        
    
    %Plot the electrodes
    fill([right_plane, right_plane, ...
        newMax, newMax], [bounds, fliplr(bounds)],[0.5, 0.5, 0.5],...
        'EdgeColor','none','FaceAlpha',0.25);
    fill([left_plane, left_plane,...
        -newMax, -newMax], [bounds, fliplr(bounds)],[0.5, 0.5, 0.5],...
        'EdgeColor','none','FaceAlpha',0.25);
    
    %plot images of image planes for reference
    for i = 1:nImageIterations
        l_image_l = jth_left_image_Xcoords(right_plane,left_plane,right_plane,i);
        l_image_r = jth_right_image_Xcoords(right_plane,left_plane,right_plane,i);
        r_image_l = jth_left_image_Xcoords(left_plane,left_plane,right_plane,i);
        r_image_r = jth_right_image_Xcoords(left_plane,left_plane,right_plane,i);
        plot([l_image_l l_image_l],bounds,'--','Color',[0 0 0]);
        plot([l_image_r l_image_r],bounds,'--','Color',[0 0 0]);
        plot([r_image_l r_image_l],bounds,'--','Color',[0 0 0]);
        plot([r_image_r r_image_r],bounds,'--','Color',[0 0 0]);
    end
    
    xlabel('X (nm)');
    ylabel('Y (nm)');
    
end