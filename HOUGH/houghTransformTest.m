%{
    Test implementation of a Hough transform
%}

function [HS] = houghTransformTest(Input,range)

    HS = 0;
    
    [sY,sX,~] = size(Input);
    
    
    [X,Y] = meshgrid(1:10,1:10);
    pixels = [X(:),Y(:)];
    
    s = size(pixels);
    for n = 1:s(1);
        [x,y] = pixels(n,:);
        
    end
end