function hotcold_heatmap_of_matrix(matrix,combineAcrossDiag,cmap)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: create a heatmap of a matrix with positive and
    %negative values by plotting filled squares with a color scale (defualt
    %is blue-to-red). If there aren't too many cells, the value of each
    %squeare will also be printed within it. 
    %
    %~~~INPUTS~~~:
    %
    %matrix: a 2D matrix of values, usually square
    %
    %combineAcrossDiag: logical variable; if true (and matrix is square),
    %   values in the (i,j) and (j,i) locations of the matrix will be added
    %   together and only plotted in the (i,j) location.
    %
    %cmap: an optional 3-column array containing a custom color map. If not
    %   included, the default blue-to-red colormap will be used. 

    
    %Default inputs
    if nargin < 3
        cmap = importdata('hot_cold_cmap.mat'); %Default color map
    end
    if nargin < 2
        combineAcrossDiag = false;
    end
   
    N1 = size(matrix,1);
    N2 = size(matrix,2);
    if combineAcrossDiag && N1 ~= N2
        combineAcrossDiag = false;
        warning('Cannot combine across diagonal with non-sqaure matrix');
    end
    
    %Combine across diagonal is requested (and possible)
    if combineAcrossDiag
        for i = 1:N1
            for j = i+1:N1
                matrix(i,j) = matrix(i,j) + matrix(j,i);
                matrix(j,i) = 0;
            end
        end
    end
    
    Nm = max(N1,N2);
    if Nm < 13
        edge_color = [0 0 0];
    else
        edge_color = 'none';
    end
      
    cbound = max(max(abs(matrix)));
    n = size(cmap,1);
    
    figure();
    hold on;
    for i = 1:N1
        for j = 1:N2
            %First let's figure out the color
            v = matrix(i,j);
            c = ceil((v + cbound)/(2*cbound)*n);
            if c > n
                c = n;
            elseif c < 1
                c = 1;
            end
            
            %Get color for square
            if combineAcrossDiag && i > j
                col = [0.5 0.5 0.5];
            else
                col = cmap(c,:);
            end
            
            %Now we can make the square
            xvals = [i i+1 i+1 i i];
            yvals = [j j j+1 j+1 j];
            fill(xvals,yvals,col,'EdgeColor',edge_color);
            
            %Don't try to plot text if there's too much of it!
            if Nm < 13 && (~combineAcrossDiag || i <= j)
                text(i+0.5,j+0.5,num2str(v),'HorizontalAlignment','center');
            end
        end
    end
    
    colormap(cmap);
    caxis([-cbound cbound]);
    colorbar();
    
    xticks((1.5:N1+0.5));
    yticks((1.5:N1+0.5));  

end