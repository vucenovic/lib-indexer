%{

%}
function [] = main(imagePath)
    image = imread(imagePath);
    [image, distortionData] = PerspectiveCorrection(image);
    
    shelf_height_px = calc_shelf_height(distortionData);
    labelQuads = label_detection(image, shelf_height_px, true);
    
    %labelQuads = [50,100,400,500;100,50,500,300;700,50,1000,300];
    
    labelQuads = sort_labels(labelQuads);
    
    labels = [];
    for i = 1:length(labelQuads)
        labels = [labels, struct("labels",[])];
        for j = 1:size(labelQuads(i).labels,1)
            label = ocr(imcrop(image,labelQuads(i).labels(j,:)));
            label.bounds = labelQuads(i).labels(j,:);
            labels(i).labels = [labels(i).labels, label];
        end
    end

    eliminate_duplicates();     % remove duplicate neighbors
    
    jsonData = jsonencode(labels);
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

%{
    Sorts the labels into rows and sorts those rows

    Author:
        Anand Eichner (11808244)
%}
function lines = sort_labels(labels)
    %% Group Labels into lines
    THRES = (labels(1,3) - labels(1,1)) * 0.75;
    lines = [struct("labels",[labels(1,:)],"average", (labels(1,1) + labels(1,3))/2)];
    for i = 2:size(labels,1)
        found = false;
        for j = 1:length(lines)
            lHeight = (labels(i,1) + labels(i,3))/2;
            if (abs(lines(j).average-lHeight) < THRES)
                lines(j).average = (lines(j).average * size(lines(j).labels,1) + lHeight) / (size(lines(j).labels,1)+1);
                lines(j).labels = [lines(j).labels; labels(i,:)];
                found = true;
                break;
            end
        end
        if (found == false)
            lines = [lines, struct("labels",[labels(i,:)],"average", (labels(i,1) + labels(i,3))/2)];
        end
    end
    
    %% Sort lines in ascending order
    for i = 1:length(lines)
        lines(i).labels = sortrows(lines(i).labels,2);
    end
end