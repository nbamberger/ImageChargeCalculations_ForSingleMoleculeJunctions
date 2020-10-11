function Xcoords = jth_left_image_Xcoords(Xcoords, left, right, iterNum)
    %Copyright 2020 LabMonti.  Written by Nathan Bamberger.  This work is 
    %licensed under the Creative Commons Attribution-NonCommercial 4.0 
    %International License. To view a copy of this license, visit 
    %http://creativecommons.org/licenses/by-nc/4.0/.  
    %
    %Function Description: Implements closed-form solution for position of
    %jth images on the left (only applied to x-coordinates)
    %
    %~~~INPUTS~~~:
    %
    %Xcoords: vector of x coordinates for a set of points between the two
    %   electrodes
    %
    %left/right: x-coordinate of left/right image plane
    %
    %iterNum: what order image charges are being considered, A.K.A j of
    %   "jth images"
    %
    %######################################################################
    %
    %~~~OUTPUTS~~~:
    %   
    %Xcoords: vector of x coordinates for the jth images on the left of the
    %   original points
    
    
    Xcoords = Xcoords*(-1)^iterNum + 2*ceil(iterNum/2)*left ...
        -2*floor(iterNum/2)*right;

end