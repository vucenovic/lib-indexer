img = imread('input2.jpg');
img_double = im2double(img);
img_grey = rgb2gray(img_double);

img_threshold = (img_grey > 0.85) & (img_grey < 0.95);

integral_img = integralImage(img_threshold);

%{
imshow(img_threshold);
figure(2);
se = strel('rectangle', [5,5]);
imshow(imdilate(img_threshold, se));
%}

%% CORNER DETECTION
img_grey_gaus = imgaussfilt(img_grey, 4);
corners = corner(img_grey_gaus, 5000, 'FilterCoefficients', fspecial('gaussian',[9 1], 5));
img_corner_highlight = img_grey;

imshow(img_corner_highlight);
hold on;
for c = 1:size(corners, 1)
    th = 0:pi/50:2*pi;
    r = 5;
    x = r * cos(th) + corners(c, 1);
    y = r * sin(th) + corners(c, 2);
    plot(x, y, 'r');
end
hold off;

%% INTEGRAL IMAGING
img_inverse = imcomplement(img_grey);

integral_img_brights = integralImage(img_grey);
integral_img_darks = integralImage(img_inverse);

%% FIND LABELS
brights_darks_ratio_range = [10, 20];   % ratio between II brightness sum/II darkness sum [min acceptable, max acceptable]
x_diff_range = [50, 100];               % absolute distance that two corners need to be apart in x direction [min, max]
y_diff_range = [75, 125];               % absolute distance that two corners need to be apart in y direction [min, max]
labels = [];

for c = 1:size(corners, 1)              % for every corner, find neighbors within x_diff_range and y_diff_range and check bright/dark ratio
    % quadtree to find neaarest neighbor with minimum distance of
    % [x_diff_range(1), y_diff_range(1)] and maximum distance of [x_diff_range(2), y_diff_range(2)]
    % if no other corner is in range, ignore the current corner and move on
    % if another corner is found in range, do integral imaging and
    % calculate the brightness/darkness ratio

    neighbors = find_neighbors(corners, corners(c), x_diff_range, y_diff_range);
    for n = 1:size(neighbors, 1)
        if check_if_label(corners(c), neighbors(1), brights_darks_ratio_range)
            labels = [labels; neighbors(1)];
        end
    end
    
end

% labels contains final result

%% FUNCTIONS
function result = find_neighbors(coords, target_coords, x_diff_range, y_diff_range)
    % find nearest neighbors in range (e.g. with quadtree)
end

function result = check_if_label(target_coords, target_coords_2, brights_darks_ratio_range)
    % integral imaging
end

