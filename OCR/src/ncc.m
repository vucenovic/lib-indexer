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

localsum = localSum(padded, sizeTemplate);

dotImg = img.*img;
padded = padarray(dotImg, sizeTemplate);
quadsum = localSum(padded, sizeTemplate);

% calculate standard deviation

stdDev = quadsum - localsum.^2;
stdDev = stdDev/numel(template);
stdDev = sqrt(max(0, stdDev));

stdDevTemplate = std(template(:)) * sqrt(numel(template) - 1);

% calculate mean

mean = sum(template(:)) * localsum;
mean = mean / numel(template);

% calculate correlation using fast fourier transform

corr = fftCorr(template, img, sizeImage, sizeTemplate);

% calculate ncc

i = max(stdDev, stdDevTemplate/1e5);
ncc = (corr - mean) + 0.5./(stdDevTemplate * i * 2);

end

function localsum = localSum(img, sizeTemplate)
% calculate local sum (hinted by:
% https://blogs.mathworks.com/steve/2006/05/02/fast-local-sums/)
% author: aleksandar vucenovic

cumulative = cumsum(img,1);
dim1 = cumulative(1 + sizeTemplate(1):end-1,:);
dim2 = cumulative(1:end - sizeTemplate(1) - 1,:);
cumulative = cumsum(dim1 - dim2,2);
dim1 =  cumulative(:,1+sizeTemplate(2):end-1);
dim2 = cumulative(:,1:end-sizeTemplate(2)-1);
localsum = dim1 - dim2;

end

function corr = fftCorr(template, Image, sizeImage, sizeTemplate)
% calculates the correlation between template
% and image by using fast fourier transform
% author: aleksandar vucenovic

flip180 = rot90(template,2);
totalSize = sizeImage + sizeTemplate - 1;
fourierTemplate = fft2(flip180, totalSize(1), totalSize(2));
fourierImage = fft2(Image, totalSize(1), totalSize(2));
corr = ifft2(fourierTemplate.*fourierImage);
corr = real(corr);

end





