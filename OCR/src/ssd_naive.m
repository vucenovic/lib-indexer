% SUM OF SQUARED DIFFERENCES ALGORITHM (NAIVE)
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% a naive implementation of the SSD algorithm, which doesnt use xcorr
% and only the sum of squared differences. this algorithm makes use 
% of the template and blob (potential character) being the same size
% and naively overlaps them and calculates the SDD. 

% the idea of SDD is, that the sum of squared differences is low (close to
% 0) if the images correlate. if you wish to use the non-naive version
% choose "ssd" instead.

% result: a minimum double value

% usage: use template and image 

% example:

% using A.bmp and "A" found in label
% ssd: 282
% using L.bmp and "A" found in label
% max correlation: 520 

function ssd = ssd_naive(template, image)
ssd = 0;
template = imcomplement(imread('temp/A.bmp'));
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

if (size(template) > size(img))
    error('template must be equal to image!')
end

if (size(template) < size(img))
    error('template must be equal to image!')
end

% calc SD for each pixel and sum up

for i = 1:size(img,1)
    for j = 1:size(img,2)
        squarediff = (img(i,j) - template(i,j)).^2;
        ssd = ssd + squarediff;
    end
end

end