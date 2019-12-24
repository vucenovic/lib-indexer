% PREPROCESSING FOR OCR 
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% Preprocessing steps for the OCR are done in this file. Preprocessing
% steps include image upscaling,for a better letter and digit quality,
% image straightening, for the template matching not having to rotate the
% templates around the image, boundingboxes and finally the image
% segmentation, to isolate each letter and digit for the template matching.

% return: an array of subImages, which are segmented digits/letters or
% lines

% usage: call preprocessing to get an array of subImages, which can be
% template matched 

% todo: resize subImages to optimize and enable template matching

function img = preprocessing(~)
% convert to binary image

img = imread('../test/Label_1.png');
[x, y, z] = size(img); 

% upscale and apply adaptive threshold

img = imresize(img, 3);

if z == 3 
   img = rgb2gray(im2double(img)); 
end

img = imbinarize(img,'adaptive','ForegroundPolarity','dark','Sensitivity',0.45);

% straighten image

angle = calcAngle(img, 0.1);
imshowpair(imrotate(img, -angle, 'bicubic'), img, 'montage');
img = imrotate(img, -angle, 'bicubic');

imshow(img);

% dilate and fill 

edgeImg = edge(img, 'prewitt');
se = strel('square',2);                 % structuring element for dilation
edgeImgDilate = imdilate(edgeImg, se); 
imshow(edgeImgDilate);
filledImg= imfill(edgeImgDilate,'holes');
imshow(filledImg);

% use regionprops to get bounding boxes of objects

box = regionprops(logical(filledImg), 'BoundingBox');
imshow(img);
hold on;
colors = hsv(numel(box));
for k = 1:numel(box)
    rectangle('position',box(k).BoundingBox, 'EdgeColor',colors(k,:));
end

% segment the regions by cropping image using bounding box rectangle
% coordinates

img = {};

for k = 3:numel(box)
    subImage = imcrop(img, box(k).BoundingBox);
    imshow(subImage);
    img(end+1) = subImage;
end
end


function [angle] = calcAngle(image, precision)
% calculates the angle to be rotated at in a range of -45 to 45 degrees
% usage: calculate the angle and rotate the image using 'bicubic' method
% author: aleksandar vucenovic, 01635282

% precondition

if ~ismatrix(image)
    error('The image must be binarized before applying this function!')
end

% angle calculation using hough

angle = angleHough(image, precision);

angle = mod(45+angle,90)-45;            
end


function angle = angleHough(image, precision)
    % edge detection using prewitt
    
    BW = edge(image,'prewitt');
    
    % perform the hough transform.
    
    [H, T, ~] = hough(BW,'Theta',-90:precision:90-precision);
    
    % find the dominant lines, by calculating variance at angles, folding
    % image, return column to angle 
    
    data=var(H);                      
    fold=floor(90/precision);         
    data=data(1:fold) + data(end-fold+1:end);
    [~, column] = max(data);          
    angle = -T(column);               
end