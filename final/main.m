%{

%}
function jsonData = main(imagePath)
    baseimage = imread(imagePath);
    %% fix jpg orientation
    info = imfinfo(imagePath);
    if isfield(info,'Format') && info(1).Format == "jpg" && isfield(info,'Orientation')
       orient = info(1).Orientation;
       switch orient
         case 1
            %normal, leave the data alone
         case 2
            baseimage = baseimage(:,end:-1:1,:);         %right to left
         case 3
            baseimage = baseimage(end:-1:1,end:-1:1,:);  %180 degree rotation
         case 4
            baseimage = baseimage(end:-1:1,:,:);         %bottom to top
         case 5
            baseimage = permute(baseimage, [2 1 3]);     %counterclockwise and upside down
         case 6
            baseimage = rot90(baseimage,3);              %undo 90 degree by rotating 270
         case 7
            baseimage = rot90(baseimage(end:-1:1,:,:));  %undo counterclockwise and left/right
         case 8
            baseimage = rot90(baseimage);                %undo 270 rotation by rotating 90
         otherwise
            warning(sprintf('unknown orientation %g ignored\n', orient));
       end
    end
    
    %% actual stuff
    tic
    [image, distortionData] = PerspectiveCorrection(baseimage);
    toc
    
    shelf_height_px = calc_shelf_height(distortionData);
    labelQuads = label_detection(image, shelf_height_px, false);
    
    labelQuads = sort_labels(labelQuads);
    toc
    
    labels = [];
    for i = 1:length(labelQuads)
        labels = [labels, struct("labels",[])];
        for j = 1:size(labelQuads(i).labels,1)
            labelQuad = labelQuads(i).labels(j,:);
            label = ocrWrapper(image(labelQuad(2):labelQuad(4), labelQuad(1):labelQuad(3)));
            label.bounds = backcorrectLabel(labelQuads(i).labels(j,:),distortionData);
            labels(i).labels = [labels(i).labels, label];
        end
    end
    toc
    
    labels = remove_invalid_labels(labels);
    
    jsonData = jsonencode(labels);
    
    displayOutput(labels,baseimage);
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
    A wrapper around the builtin ocr to convert its results to our format
    and execute sanity checks

    Takes: an image of the label
%}
function label = ocrWrapper(subImage)
    result = ocr(subImage);
    stringsRaw = splitlines(string(result.Text));
    strings = [];
    word = 1;
    for i = 1:length(stringsRaw)
        str = stringsRaw(i);
        if word == 1
            r = regexp(str, "[A-Za-z]{3}");
            if length(r) == 1
                strings = [strings,str];
                word = 2;
            end
        elseif word == 2
            r = regexp(str, "[\dO]{3}");
            if length(r) == 1
                strings = [strings,str];
                word = 3;
            end
        else
            r = regexp(str, "[A-Za-z]{3,}");
            if length(r) == 1
                strings = [strings,str];
                word = 4;
            end
        end
    end
    
    if length(strings)<3
        label.wordOne = "";
        label.wordTwo = "";
        label.author = "";
    else
        label.wordOne = strings(1);
        label.wordTwo = strings(2);
        label.author = strings(3);
    end
end

%{
    Takes the label bounds outputted by the label detection and
    backcorrects them to a bounding quad in the original image

    Takes:
%}
function quad = backcorrectLabel(labelQuad,spaceData)
    points = [labelQuad([1,2]);labelQuad([3,2]);labelQuad([3,4]);labelQuad([1,4])];
    imOffset = [spaceData.imRef.XWorldLimits(1),spaceData.imRef.YWorldLimits(1)];
    points = points + repmat(imOffset,4,1);
    quad = transformPointsInverse(spaceData.transform,points);
    quad = round(quad);
end

%{
    Visualizes our results
%}
function [] = displayOutput(labels,image)
    figure;
    imshow(image), hold on;
    for i = 1:length(labels)
        label = labels(i);
        clr = sscanf("6DBDC9",'%2x%2x%2x',[1 3])/255;
        plot(polyshape(label.bounds(:,1),label.bounds(:,2)),"FaceAlpha",1, "FaceColor", "white", "EdgeColor", clr, "LineWidth", 2);
        t = text((label.bounds(1,1) + label.bounds(3,1))/2,(label.bounds(1,2) + label.bounds(3,2))/2,...
            [label.wordOne,label.wordTwo,label.author],...
            "HorizontalAlignment","center",'FontSize',10);
        uistack(t, 'top');
    end
end

%{
    Book labels only have at max 8 chars for their author.
    So we eiminate all labels with more chars than that. They are no valid
    labels.

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = remove_invalid_labels(labels)
    
    result = [];
    for i = 1:size(labels, 2)
        for j = 1:size(labels(i).labels, 2)
            label = labels(i).labels(j);
            if strlength(label.author) <= 8 && ...
               strlength(label.wordOne) == 3 && ...
               strlength(label.wordTwo) == 3 && ...
               strlength(label.author) > 0
                result = [result; label];
            end
        end
    end

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