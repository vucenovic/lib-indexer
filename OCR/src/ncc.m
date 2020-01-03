% NORMALIZED CORRELATION COEFFICIENT ALGORITHM
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% Optical character recognition through template matching using
% NCC-Algorithm. The normalized correleation coefficient algorithm
% determines the corelation of pixels between two images and returns a
% value between -1 and 1, wheras 1 means they fully correlate, -1 means it
% correlates with its negation, 0 means no correlation.

% NCC description:
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

% calculate local sum (hinted by:
% https://blogs.mathworks.com/steve/2006/05/02/fast-local-sums/)

localsum = localSum(img, sizeTemplate);

quadImg = img.^2;
quadsum = localSum(quadImg, sizeTemplate);

% calculate standard deviation and mean

diff = quadsum - localsum.^2/numel(template);
sigmaImg = sqrt(max(0, diff));
stdTemplate = std(template(:));
meanTemplate = mean(mean(template));

if stdTemplate == 0
    error("the values of the template must not be the same!")
end

sigmaTemplate = sqrt(numel(template) - 1) * stdTemplate;
denom = sigmaTemplate * sigmaImg;

% calculate correlation using fast fourier transform

corr = fftCorr(template, img, sizeImage, sizeTemplate);
numer = (corr - localsum * meanTemplate);

% calculate ncorr

correlation = numer/denom;

% set values outside of -1 and 1 to zero (happens through variance)


end

function localsum = localSum(img, sizeTemplate)
% calculates the correlation between template
% and image by using fast fourier transform
% author:  Dirk-Jan Kroon (2020). Fast/Robust Template Matching
%(https://www.mathworks.com/matlabcentral/fileexchange/24925-fast-robust-template-matching),
% MATLAB Central File Exchange. Retrieved January 3, 2020. 

% note: minor changes to author's code for readability

cumulative = cumsum(img,1);
dim1 = cumulative(1 + sizeTemplate(1):end-1,:);
dim2 = cumulative(1:end - sizeTemplate(1) - 1,:);
cumulative = cumsum(dim1 - dim2,2);
dim1 =  cumulative(:,1+sizeTemplate(2):end-1);
dim2 = cumulative(:,1:end-sizeTemplate(2)-1);
localsum = dim1 - dim2;

end

function corr = fftCorr(template, Image, sizeImage, sizeTemplate)
% author: aleksandar vucenovic

% note: followed this paper 
% http://scribblethink.org/Work/nvisionInterface/nip.pdf

flip180 = rot90(template,2);
totalSize = sizeImage + sizeTemplate - 1;
fourierTemplate = fft(flip180, totalSize(1), totalSize(2));
fourierImage = fft(Image, totalSize(1), totalSize(2));
corr = ifft(fourierTemplate.*fourierImage);
corr = real(corr);

end

function result = forceSize(mat, target)
one = target;

for k = 1:2
    if size(mat,k) > target[k]
        difference = (size(mat,k) - target[k])/2;
        one = [floor:
end
end





