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

% return: a correlation coefficient matrix, which values range from -1 to 1

% usage: input a template and image 

% examples:

% using A.bmp and "A" found in label
% max correlation: 0.5324 
% using L.bmp and "A" found in label
% max correlation: 0.3155

function correlation = ncc(template, img)
template = imcomplement(template);
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

correlation = numer./denom;

% set values outside of -1 and 1 to zero (happens through variance)

correlation(correlation > 1) = 0;
correlation(correlation < -1) = 0;

end

function localsum = localSum(I,T_size)
% author:  Dirk-Jan Kroon (2020). Fast/Robust Template Matching
%(https://www.mathworks.com/matlabcentral/fileexchange/24925-fast-robust-template-matching),
% MATLAB Central File Exchange. Retrieved January 3, 2020. 

B = padarray(I,T_size);
% Calculate for each pixel the sum of the region around it,
% with the region the size of the template.
if(length(T_size)==2)
    % 2D localsum
    s = cumsum(B,1);
    c = s(1+T_size(1):end-1,:)-s(1:end-T_size(1)-1,:);
    s = cumsum(c,2);
    localsum = s(:,1+T_size(2):end-1)-s(:,1:end-T_size(2)-1);
end

end

function corr = fftCorr(template, Image, sizeImage, sizeTemplate)
% author: aleksandar vucenovic

% note: followed this paper 
% http://scribblethink.org/Work/nvisionInterface/nip.pdf

flip180 = rot90(template,2);
totalSize = sizeImage + sizeTemplate - 1;
fourierTemplate = fft2(flip180, totalSize(1), totalSize(2));
fourierImage = fft2(Image, totalSize(1), totalSize(2));
corr = ifft2(fourierTemplate.*fourierImage);
corr = real(corr);
corr = corr(1:totalSize(1),1:totalSize(2));

end







