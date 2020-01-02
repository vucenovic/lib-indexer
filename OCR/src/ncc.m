% NORMALIZED CORRELATION COEFFICIENT ALGORITHM
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% Optical character recognition through template matching using
% NCC-Algorithm. The normalized correleation coefficient algorithm
% determines the corelation of pixels between two images and returns a
% value between -1 and 1, wheras 1 means they fully correlate, -1 means it
% correlates with its negation, 0 means no correlation.

% This code follows the following paper:
% http://scribblethink.org/Work/nvisionInterface/nip.pdf

% The correlation coefficient that must be reached at its peak is 0.5,
% determined by trying out different templates and images, a correct result
% was achieved when the factor was >0.5 and incorrect when <0.5

% return: a correlation coefficient, which ranges from -1 to 1

% usage: input a template and image 

function coefficient = ncc(template, img)
template = imcomplement(imread('temp/L.bmp'));
img = imread('temp/tempSubImage11.png');

% check prerequisites

[x, y, z] = size(img); 
[u, v, w] = size(template); 

if z == 3 
   img = rgb2gray(im2double(img)); 
end

if w == 3 
   img = rgb2gray(im2double(template)); 
end

if size(template) > size(img)
    error('template must be smaller or equal size of image!')
end

% compute correlations matrix with ncc

correlation = [];

sizeTemplate = size(template);
sizeImage = size(img);

% add padding to original image

padded = padarray(img, sizeTemplate);

% calculate local sum (hinted by:
% https://blogs.mathworks.com/steve/2006/05/02/fast-local-sums/)

cumulative = cumsum(padded,1);
dim1 = cumulative(1 + sizeTemplate(1):end-1,:);
dim2 = cumulative(1:end - sizeTemplate(1) - 1,:);
cumulative = cumsum(dim1 - dim2,2);
dim1 =  cumulative(:,1+sizeTemplate(2):end-1);
dim2 = cumulative(:,1:end-sizeTemplate(2)-1);
localsum = dim1 - dim2;

% calculate standard deviation

stdDev = 




end





