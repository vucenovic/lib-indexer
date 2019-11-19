% PREPROCESSING FOR OCR 
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% convert to binary image

img = imread('../test/Label_1.png');
[x, y, z] = size(img); 

% upscale

img = imresize(img, 3);

if z == 3 
   img = rgb2gray(im2double(img)); 
end

img = imbinarize(img,'adaptive','ForegroundPolarity','dark','Sensitivity',0.45);

% straighten image

angle = horizon(img, 0.1, 'hough');
imshowpair(imrotate(img, -angle, 'bicubic'), img, 'montage');
img = imrotate(img, -angle, 'bicubic');

% remove lines

img = edge(img, 'canny');
img = removeLines(removeLines(removeLines(removeLines(removeLines(img)))));
img = removeLines(removeLines(removeLines(removeLines(removeLines(img)))));
img = removeLines(removeLines(removeLines(removeLines(removeLines(img)))));
img = removeLines(removeLines(removeLines(removeLines(removeLines(img)))));
img = removeLines(img);

imshow(img);

function [img] = removeLines(img)
% removes horitontal lines of image
% author: aleksandar vucenovic, 01635282

[H,theta,rho] = hough(img);
P = houghpeaks(H,10,'NHoodSize',[1 1]);
lines_found = houghlines(img,theta,rho,P,...
    'FillGap',500,'MinLength',1);
for k = 1:length(lines_found)
   % extract one line:
   xy = [lines_found(k).point1; lines_found(k).point2];
   % remove the lines from the image:
   % note that I take a buffer of 3 to the 'width' of the line
   img(xy(1,2):xy(1,2)+3,xy(1,1):xy(2,1)) = 0;
   img(xy(1,2):xy(1,2)+5,xy(1,1):xy(2,1)) = 0;
end
end




function [angle] = horizon(image, varargin)
% HORIZON estimates the horizon rotation in the image.
%   ANGLE=HORIZON(I) returns rotation of an estimated horizon
%   in the image I. The returned value ANGLE is in the
%   range <-45,45> degrees.
%
%   ANGLE=HORIZON(I, PRECISION) aligns the image I with
%   the predefined precision. The default value is 1 degree. If higher
%   precision is required, 0.1 could be a good value.
%
%   Example
%   -------
%       image = imread('board.tif');
%       angle = horizon(rgb2gray(image), 0.1, 'hough')
%       imshow(imrotate(image, -angle, 'bicubic'));
%
%   The example aligns the default image in Image Processing Toolbox.
% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 3                   % only want 3 optional inputs at most
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
optargs = {1, 'fft', 2};            % set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[precision, method, diskSize] = optargs{:};  % use memorable variable names
% Check image dimension.
if ndims(image)~=2
    error('The image must be two-dimensional (i.e. grayscale).')
end
if strcmpi(method, 'hough')
    angle = horizonHough(image, precision);
end
% Return the angle
angle = mod(45+angle,90)-45;            % rotation in -45..45 range
end
function angle = horizonHough(image, precision)
    % Detect edges.
    BW = edge(image,'prewitt');
    % Perform the Hough transform.
    [H, T, ~] = hough(BW,'Theta',-90:precision:90-precision);  
    % Find the most dominant line direction.
    data=var(H);                      % measure variance at each angle 
    fold=floor(90/precision);         % assume right angles & fold data
    data=data(1:fold) + data(end-fold+1:end);
    [~, column] = max(data);          % the column with the crispiest peaks
    angle = -T(column);               % column to degrees 
end