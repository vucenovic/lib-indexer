% PREPROCESSING FOR OCR 
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% Preprocessing steps for the OCR are done in this file. Preprocessing
% steps include image upscaling,for a better letter and digit quality,
% image straightening, for the template matching not having to rotate the
% templates around the image, boundingboxes and finally the image
% segmentation, to isolate each letter and digit for the template matching.

% return: a cell array containing all the blobs found in the image

% usage: call preprocessing to get a cellarray containg all image blobs,
% which are potential characters 

function patches = preprocessing(img)
% convert to binary image

img = img;
img = imread('Label_1.png');
[x, y, z] = size(img); 

% upscale and apply adaptive threshold

img = imresize(img, 3);

if z == 3 
   img = rgb2gray(im2double(img)); 
end

img = imbinarize(img,'adaptive','ForegroundPolarity','dark','Sensitivity',0.45);

% straighten image

angle = calcRotationAngle(img);
img = imrotate(img, -angle, 'bicubic');
img = 1-imclearborder(1 - img);

% dilate and fill 

edgeImg = edge(img, 'prewitt');
se = strel('square',2);                
edgeImgDilate = imdilate(edgeImg, se); 
filledImg= imfill(edgeImgDilate,'holes');

% use regionprops to get bounding boxes of objects

box = regionprops(logical(filledImg), 'BoundingBox', 'Centroid');

% delete lines to get characters only

box_corrected = deleteLines(box);

% slice boxes which include two characters

box_sliced = sliceBoxes(box_corrected);

% get centroids of characters
centroidsXY = vertcat(box_sliced.Centroid);

% sort the indices column-wise
% source: 
% https://stackoverflow.com/questions/43076798/how-to-control-the-order-of-detected-objects-by-regionprops-in-matlab
box_sliced = sortIndex(box_sliced, centroidsXY);
[~, ~, centroidsXY(:, 2)] = histcounts(centroidsXY(:, 2), 3); 
[~, sortIndex] = sortrows(centroidsXY, [2 1]);  
box_sliced = box_sliced(sortIndex);  

imshow(img);
hold on;
colors = hsv(numel(box_sliced));
for k = 1:length(box_sliced)
    rectangle('position',box_sliced(k).BoundingBox, 'EdgeColor',colors(k,:));
end

% segment the regions by cropping image using bounding box rectangle
% coordinates, save the first three letters of a label first,
% then the 3 digit code, then the author
patches = []
for k = 1:length(box_sliced)
        subImage = imcrop(img, box_sliced(k).BoundingBox);
        patches = [patches, struct("image",subImage)];
end

end



function angle = calcRotationAngle(image)
% calculates the angle to be rotated at in a range of -45 to 45 degrees
% usage: calculate the angle and rotate the image using 'bicubic' method
% author: aleksandar vucenovic, 01635282

% precondition

if ~ismatrix(image)
    error('The image must be binarized!')
end

% angle calculation using hough

% edge detection using prewitt
    
BW = edge(image,'prewitt');
    
% perform the hough transform.
    
[H, T, ~] = hough(BW,'Theta',-90:0.1:89.9);
    
% find the dominant lines, by calculating variance at angles, folding
% image, return column to angle 
    
data = var(H);                              
data = data(1:900) + data(end-900+1:end);
[~, column] = max(data);          
angle = -T(column);             

angle = mod(45 + angle,90) - 45;            
end

function box_sliced = sliceBoxes(box_corrected)
% author: anand eichner & aleksandar vucenovic
% slice boxes which include two characters
% the characters are monospaced, and the width should always be smaller
% than the height

box_sliced = []
for i = 1:length(box_corrected)
    b = box_corrected(i);
    coord = b.BoundingBox;
    if coord(3) > coord(4)
        wide = coord(3) / 2;
        x = coord(1);
        b1.BoundingBox = [coord(1), coord(2), wide, coord(4)];
        b2.BoundingBox = [coord(1) + wide, coord(2), wide, coord(4)];
        b1.Centroid = b.Centroid;
        b1.Centroid(1) = coord(1) + wide / 2;
        b2.Centroid = b.Centroid;
        b2.Centroid(1) = coord(1) + wide / 2 * 3;
        box_sliced = [box_sliced; b1 ; b2];
    else
        box_sliced = [box_sliced; b];
    end
end
end

function box_corrected = deleteLines(box)
% author: aleksandar vucenovic
% delete lines from label to only have boxes and centroids around chars
box_corrected = [];
for i = 1:length(box)
    b = box(i);
    coord = b.BoundingBox(3:4);
    if coord(1) < (coord(2) * 10)
        box_corrected = [box_corrected; b];
    end
end
end
