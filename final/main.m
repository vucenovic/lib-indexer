%{

%}
function [] = main(imagePath)
    image = imread(imagePath);
    [image, distortionData] = PerspectiveCorrection(image);
    
    shelf_height_px = calc_shelf_height(distortionData);
    %labelQuads = label_detection(image, shelf_height_px, true);
    
    labelQuads = [50,100,700,500;100,50,1000,300];
    
    STRATEGY = "NCC";
    for i = 1:size(labelQuads, 1)
        quad = labelQuads(i, :);
        subImg = image(quad(2):quad(4),quad(1):quad(3), :);
        patches = preprocessing(subImg);
        for j = 1:length(patches)
            if(STRATEGY=="NCC")
               ncc(); 
            elseif(STRATEGY=="SSD")
                
            end
        end
    end
    
    sort_labels();              % labels to shelf coords
    eliminate_duplicates();     % remove duplicate neighbors
    
end

%{
    Find the height in px for a shelf.
    Basically takes the difference in px of the two perspective correction-
    horizontals that are furthest apart.

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function shelf_height = calc_shelf_height(distortion_data)
    
    horizontals_positions = [distortion_data.horizontals.P];
    horizontals_positions = sort(horizontals_positions(2:2:end));
    horizontals_distances = diff(horizontals_positions);
    shelf_height = max(horizontals_distances);

end
