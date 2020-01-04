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
%imshowpair(imrotate(img, -angle, 'bicubic'), img, 'montage');
img = imrotate(img, -angle, 'bicubic');
img = 1-imclearborder(1 - img);

%imshow(img);
%imwrite(img, 'temp/label.png');

% dilate and fill 

edgeImg = edge(img, 'prewitt');
se = strel('square',2);                 % structuring element for dilation
edgeImgDilate = imdilate(edgeImg, se); 
%imshow(edgeImgDilate);
filledImg= imfill(edgeImgDilate,'holes');
imshow(filledImg);

% use regionprops to get bounding boxes of objects

box = regionprops(logical(filledImg), 'BoundingBox', 'Centroid');
box_corrected = [];
for i = 1:length(box)
    b = box(i);
    coord = b.BoundingBox(3:4);
    if coord(1) < (coord(2) * 10)
        box_corrected = [box_corrected; b];
    end
end


centroidsXY = vertcat(box_corrected.Centroid);
imshow(img)
hold on
plot(centroidsXY(:,1),centroidsXY(:,2),'b*')
hold off
%imshow(img);
hold on;
colors = hsv(numel(box));
for k = 1:numel(box_corrected)
    rectangle('position',box_corrected(k).BoundingBox, 'EdgeColor',colors(k,:));
end

% segment the regions by cropping image using bounding box rectangle
% coordinates, save them as images in a temporary folder
patches = []
for k = 1:numel(box_corrected)
    subImage = imcrop(img, box_corrected(k).BoundingBox);
    %imshow(subImage);
    %filename = sprintf('temp/tempSubImage%d.png', k);
    %imwrite(imresize(subImage, [42, 24]), filename);
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