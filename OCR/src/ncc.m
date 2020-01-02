% NORMALIZED CORRELATION COEFFICIENT ALGORITHM
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% Optical character recognition through template matching using
% NCC-Algorithm. The normalized correleation coefficient algorithm
% determines the corelation of pixels between two images and returns a
% value between -1 and 1, wheras 1 means they fully correlate, -1 means it
% correlates with its negation, 0 means no correlation.

% return: a correlation coefficient, which ranges from -1 to 1

% usage: input a template and image 

template = imcomplement(imread('temp/A.bmp'));
img = imread('temp/tempSubImage9.png');
c = normxcorr2(imcomplement(imread('temp/A.bmp')), imread('temp/tempSubImage9.png'));

[ypeak, xpeak] = find(c==max(c(:)));

yoffSet = ypeak-size(template,1);
xoffSet = xpeak-size(template,2);

figure
imshow(img);
imrect(gca, [xoffSet+1, yoffSet+1, size(template,2), size(template,1)]);





