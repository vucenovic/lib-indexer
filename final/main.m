%{

%}
function [] = main(imagePath)
    image = imread(imagePath);
    %% fix jpg orientation
    info = imfinfo(imagePath);
    if isfield(info,'Format') && info(1).Format == "jpg" && isfield(info,'Orientation')
       orient = info(1).Orientation;
       switch orient
         case 1
            %normal, leave the data alone
         case 2
            image = image(:,end:-1:1,:);         %right to left
         case 3
            image = image(end:-1:1,end:-1:1,:);  %180 degree rotation
         case 4
            image = image(end:-1:1,:,:);         %bottom to top
         case 5
            image = permute(image, [2 1 3]);     %counterclockwise and upside down
         case 6
            image = rot90(image,3);              %undo 90 degree by rotating 270
         case 7
            image = rot90(image(end:-1:1,:,:));  %undo counterclockwise and left/right
         case 8
            image = rot90(image);                %undo 270 rotation by rotating 90
         otherwise
            warning(sprintf('unknown orientation %g ignored\n', orient));
       end
    end
    
    %% actual stuff
    [image, distortionData] = PerspectiveCorrection(image);
    
    shelf_height_px = calc_shelf_height(distortionData);
    labelQuads = label_detection(image, shelf_height_px, true);
    
    %labelQuads = [50,100,400,500;100,50,500,300;700,50,1000,300];
    
    labelQuads = sort_labels(labelQuads);
    
    labels = [];
    for i = 1:length(labelQuads)
        labels = [labels, struct("labels",[])];
        for j = 1:size(labelQuads(i).labels,1)
            labelQuad = labelQuads(i).labels(j,:);
            label = ocr(image(labelQuad(2):labelQuad(4), labelQuad(1):labelQuad(3)));
            label.bounds = labelQuads(i).labels(j,:);
            labels(i).labels = [labels(i).labels, label];
        end
    end
    
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
    THRES = (labels(1,4) - labels(1,2)) * 0.75;
    lines = [struct("labels",[labels(1,:)],"average", (labels(1,2) + labels(1,4))/2)];
    for i = 2:size(labels,1)
        found = false;
        for j = 1:length(lines)
            lHeight = (labels(i,2) + labels(i,4))/2;
            if (abs(lines(j).average-lHeight) < THRES)
                lines(j).average = (lines(j).average * size(lines(j).labels,1) + lHeight) / (size(lines(j).labels,1)+1);
                lines(j).labels = [lines(j).labels; labels(i,:)];
                found = true;
                break;
            end
        end
        if (found == false)
            lines = [lines, struct("labels",[labels(i,:)],"average", (labels(i,2) + labels(i,4))/2)];
        end
    end
    
    %% Sort lines in ascending order
    for i = 1:length(lines)
        lines(i).labels = sortrows(lines(i).labels,1);
    end
end