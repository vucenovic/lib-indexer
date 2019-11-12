img = imread('input2.jpg');
img_double = im2double(img);
img_grey = rgb2gray(img_double);

test_img_threshold = (img_grey > 0.85) & (img_grey < 0.95);

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
img_integral = generate_integral_image(img_grey);

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
        if check_if_label(img_integral, corners(c), neighbors(1), brights_darks_ratio_range)
            labels = [labels; corners(c), neighbors(1)];
        end
    end
    
end

% labels contains final result

%% FUNCTIONS
function result = find_neighbors(coords, target_coords, x_diff_range, y_diff_range)
    % find nearest neighbors in range (e.g. with quadtree)
end

%{
    Returns true if the two target vectors could describe a label (based on
    a heuristic), or false if it is no label that can be detected.

    Sources:
        http://www.florian-oeser.de/wordpress/wp-content/2012/10/crow-1984.pdf
        accessed on 2019/11/12
    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = check_if_label(integral_image, target_coords, target_coords_2, brights_darks_ratio_range)
    
    A = integral_image(target_coords);
    B = integral_image([target_coords(1), target_coords_2(2)]);
    C = integral_image([target_coords(2), target_coords_2(1)]);
    D = integral_image(target_coords_2);
    
    target_diff   = diff([target_coords(1), target_coords_2(1), target_coords(2), target_coords_2(2)]);
    ii_max_value  = target_diff(1) * target_diff(2);
    ii_brightness = D + A - B - C;
    ii_darkness   = ii_max_value - ii_brightness;
    
    ii_bright_dark_ratio = ii_brightness / ii_darkness;
    
    result = ii_bright_dark_ratio >= brights_darks_ratio_range(1) &
             ii_bright_dark_ratio <= brights_darks_ratio_range(2);

end

%{
    Returns an integral image (also called summed area table) from the
    given grayscale image (2-dimensional, [0, 1]-image).

    Sources:
        http://www.florian-oeser.de/wordpress/wp-content/2012/10/crow-1984.pdf
        accessed on 2019/11/12
    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = generate_integral_image(greyscale_image)
    
    result = greyscale_image;
    for x = 1:size(greyscale_image, 1)
       for y = 1:size(greyscale_image, 2)
           
           result(x, y) =   result( max(x-1, 1),    max(y-1, 1)) + ...
                            result(           x,    max(y-1, 1)) + ...
                            result( max(x-1, 1),              y) + ...
                            result(           x,              y);
           
       end
    end

end

