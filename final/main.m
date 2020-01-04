%{

%}
function [] = main(imagePath)
    image = imread(imagePath);
    [image,distortionData] = PerspectiveCorrection(image);
    %labelQuads = label_detection(image,false);
    
    labelQuads = [50,100,700,500;100,50,1000,300];
    
    STRATEGY = "NCC";
    for i = 1:size(labelQuads,1)
        quad = labelQuads(i,:);
        subImg = image(quad(2):quad(4),quad(1):quad(3),:);
        patches = preprocessing(subImg);
        for j = 1:length(patches)
            if(STRATEGY=="NCC")
               ncc(); 
            elseif(STRATEGY=="SSD")
                
            end
        end
    end
end