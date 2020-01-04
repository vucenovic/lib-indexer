img = imread('input2.jpg');
img_double = im2double(img);
img_grey = rgb2gray(img_double);

%% GLOBAL THRESHOLD
% we apply the threshold twice. all pixels lower than the first th are
% clipped to that same th value. this makes the second th less likely to
% be too high, yet it is high enough to separate text from labels.
 
th_img = img_grey;
global_th = otsu_threshold(th_img);
th_img(th_img < global_th) = global_th;
global_th = otsu_threshold(th_img);
th_mask = img_grey >= global_th;

%th_img = transform_clip_image(img_grey, global_th, 1);


%% CORNER DETECTION

img_grey_gaus = imgaussfilt(img_grey, 4);
corners = corner(img_grey_gaus, 5000, 'FilterCoefficients', fspecial('gaussian',[9 1], 2));
corners = [corners(:, 2), corners(:, 1)]; % swap columns so it's right
%corners = harris_corners(img_grey_gaus);


%% INTEGRAL IMAGING
% create integral image for brightness values.
% pixels higher than the global th are assumed to be white

bright_th = img_grey;
bright_th(bright_th >= global_th) = 1;
img_integral_brights = generate_integral_image(bright_th);


%% FIND LABELS

brightness_amount_range = [7, 25];      % acceptable ratio between bright and dark areas of a label (see check_if_label())
x_diff_range = [20, 125];               % distance in px that two corners need to be apart in x direction [min, max]
y_diff_range = [100, 200];              % distance in px that two corners need to be apart in y direction [min, max]
labels = [];
corners = sortrows(corners);            % sort corners by their y-axis

for c = 1:size(corners, 1)
    
    neighbors = find_neighbors(corners, c, x_diff_range, y_diff_range);
    for n = 1:size(neighbors, 1)
        if check_if_label(img_grey, th_mask, img_integral_brights, corners(c, :), neighbors(n, :), brightness_amount_range)
            labels = [labels; corners(c, :), neighbors(n, :)];
        end
    end
    
end

% labels contains final result


%% DEBUG

corners_t = [corners(:, 2), corners(:, 1)];
imshow(img_grey);
hold on;
for t = 1:size(labels, 1)
    plot([labels(t, 2), labels(t, 4)], [labels(t, 1), labels(t, 3)], 'LineWidth', 2);
end
for c = 1:size(corners_t, 1)
    th = 0:pi/2:2*pi;
    r = 5;
    x = r * cos(th) + corners_t(c, 1);
    y = r * sin(th) + corners_t(c, 2);
    plot(x, y, 'r');
end
hold off;


%% FUNCTIONS

%{

### TODO REMOVE ###

    Takes the input image (img) and transforms it so low is 0 and high is
    1.
    All values < 0 or > 1 are clipped.

    Sources:
        EVC UE Ex.3 resources: evc_histogram_clipping.m

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = transform_clip_image(img, low, high)
    result = (img - low) ./ (high - low);
    result(result < 0) = 0;
    result(result > 1) = 1;
end

%{
    Find neighboring corners that satisfy our distance diff range.
    We only look for corners that occur after coords in the coords array,
    since we already checked all neighbors that came before (we assume
    coords is sorted by y-axis).
    The resulting neighbors might form a label together with coords.

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = find_neighbors(coords, coords_target_index, x_diff_range, y_diff_range)
    
    result = [];
    for j = coords_target_index + 1:size(coords, 1)
        neighbor = coords(j, :);
        if check_for_range(coords(coords_target_index, :), neighbor, x_diff_range, y_diff_range)
            result = [result; neighbor];
        end
    end
    
end

%{
    Check if the two points (coordA, coordB) satisfy x_diff_range and
    y_diff_range.
    Labels are expected to have a certain width & height, so we restrict
    the possible corner-combinations.

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = check_for_range(coordA, coordB, x_diff_range, y_diff_range)

    diff_x = abs(coordA(2) - coordB(2));
    diff_y = abs(coordA(1) - coordB(1));
    
    result = diff_x >= x_diff_range(1) && ...
             diff_x <= x_diff_range(2) && ...
             diff_y >= y_diff_range(1) && ...
             diff_y <= y_diff_range(2);
   
end

%{
    Returns true if the two target vectors could describe a label (based on
    a heuristic), or false if it is no label that can be detected.

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = check_if_label(img, img_mask, integral_image_brights, corner, neighbor, brights_darks_ratio_range)
    
    % ### TODO REMOVE ### label_candidate = img(min([corner(1), neighbor(1)]):max([corner(1), neighbor(1)]), min([corner(2), neighbor(2)]):max([corner(2), neighbor(2)]));
    
    label_candidate_th = img_mask(min([corner(1), neighbor(1)]):max([corner(1), neighbor(1)]), min([corner(2), neighbor(2)]):max([corner(2), neighbor(2)]));
    dark_amount = sum(1-clear_label_border(label_candidate_th), 1:2);

    ii_brightness_amount = integral_image_result(integral_image_brights, corner, neighbor);    
    brightness_ratio = ii_brightness_amount / dark_amount;
    
    result = brightness_ratio >= brights_darks_ratio_range(1) && ...
             brightness_ratio <= brights_darks_ratio_range(2);

end

%{
    Clear dark areas (connected components) that touch the images' border using
    MATLAB's imclearborder.

    If more than 25% of the label is changed, we ignore the result and hand
    back the original input.
    If the center 50% (50% x, 50% y) is changed at all, we ignore the
    result and hand back the original input as well.

    Sources:
        https://blogs.mathworks.com/steve/2007/09/04/clearing-border-components/
        accessed on 2020/01/04

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = clear_label_border(label_binary)

    label_binary_inverted = 1-label_binary; % before: bright areas = 1; now: dark areas = 1
    label_binary_inverted_cleared = imclearborder(label_binary_inverted);
    label_binary_cleared = 1-label_binary_inverted_cleared;
    
    [dimensions_y, dimensions_x] = size(label_binary);
    pixel_amount = dimensions_x * dimensions_y;
    label_binary_center = label_binary(round(dimensions_x*0.25):round(dimensions_x*0.75));
    label_binary_cleared_center = label_binary_cleared(round(dimensions_x*0.25):round(dimensions_x*0.75));
    
    pixels_changed = abs(sum(label_binary, 1:2) - sum(label_binary_cleared, 1:2)); % TODO maybe swap with II
    pixels_changed_center = abs(sum(label_binary_center, 1:2) - sum(label_binary_cleared_center, 1:2));
    
    if pixels_changed > pixel_amount * 0.25 || pixels_changed_center > 0
        result = label_binary;
    else
        result = label_binary_cleared;
    end
    
end

%{
    Calculate and return the sum of intensity levels contained within the two
    given corners (must be diagonal corners; either top-left:bottom-right
    (and vice-versa) or top-right:bottom-left (and vice-versa).

    Sources:
        http://delivery.acm.org/10.1145/810000/808600/p207-crow.pdf
        accessed on 2019/11/12

        https://en.wikipedia.org/wiki/Summed-area_table#/media/File:Summed_area_table.png
        used only for variable naming
        accessed on 2019/11/12

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = integral_image_result(integral_image, corner, neighbor)

    ii_coords_1 = integral_image(corner(1), corner(2));
    ii_coords_12 = integral_image(corner(1), neighbor(2));
    ii_coords_2 = integral_image(neighbor(1), neighbor(2));
    ii_coords_21 = integral_image(neighbor(1), corner(2));
    
    % neighbor lies to the bottom left
    if corner(1) < neighbor(1) && corner(2) > neighbor(2)
        A = ii_coords_12;
        B = ii_coords_1;
        C = ii_coords_2;
        D = ii_coords_21;
        
    % neighbor lies to the bottom right
    elseif corner(1) < neighbor(1) && corner(2) < neighbor(2)
        A = ii_coords_1;
        B = ii_coords_12;
        C = ii_coords_21;
        D = ii_coords_2;
        
    else
        error('We assume the corners are sorted by Y-axis and neighbors have a Y greater than and X not equal to corner.');
    end
    
    result = D + A - B - C;
    
end

%{
    Returns an integral image (also called summed area table) from the
    given grayscale image (2-dimensional, [0, 1]-image).

    Sources:
        http://delivery.acm.org/10.1145/810000/808600/p207-crow.pdf
        accessed on 2019/11/12

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = generate_integral_image(greyscale_image)
    
    result = greyscale_image;
    for x = 1:size(greyscale_image, 1)
       for y = 1:size(greyscale_image, 2)
           
           r = 0;
           if x > 1
               r = r + result(x-1, y);
           end
           
           if y > 1 
               r = r + result(x, y-1);
           end
           
           if x > 1 && y > 1
               r = r - result(x-1, y-1);
           end
           
           result(x, y) = r + greyscale_image(x, y);
           
       end
    end

end

%{
    Calculates the greyscale_img's threshold using otsu's method.

    Sources:
        https://engineering.purdue.edu/kak/computervision/ECE661.08/OTSU_paper.pdf
        accessed on 2020/01/04

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = otsu_threshold(greyscale_img)
    
    [bin_amounts, bins] = imhist(greyscale_img);
    
    bin_size = size(bins, 1);

    % Make counts a double column vector
    bin_amounts = double(bin_amounts);

    % Variables names are chosen to be similar to the formulas in
    % the Otsu paper.
    p = bin_amounts / sum(bin_amounts);
    omega = cumsum(p);
    mu = cumsum(p .* (1:bin_size)');
    mu_t = mu(end);

    sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

    % Find the location of the maximum value of sigma_b_squared.
    % The maximum may extend over several bins, so average together the
    % locations.  If maxval is NaN, meaning that sigma_b_squared is all NaN,
    % then return 0.
    maxval = max(sigma_b_squared);
    isfinite_maxval = isfinite(maxval);
    if isfinite_maxval
        idx = mean(find(sigma_b_squared == maxval));
        % Normalize the threshold to the range [0, 1].
        t = (idx - 1) / (bin_size - 1);
    else
        t = 0.0;
    end

result = t;
% TODO ### rewrite ###

end


%{
Harris corner detector. Returns corners.
Author: Anna Berezhinskaya
Quelle: https://en.wikipedia.org/wiki/Harris_Corner_Detector
%}
function result = harris_corners(greyscale_image)
    % 2nd Step:
    %{    
     The sum of squared differences between 2 patches will be given 
     with following equation:
     (1) F(u,v) = sum([I(x+u,y+v) - I(x,y)]^2).
     Now using Taylor expansion we can approximate I(x+u,y+v) as following:
     I(x+u,y+v) = I(x,y) + dx * d(I(x,y))/dx + dy * d(I(x,y))/dy

     We can than rewrite the equation (1)  in the following way :
        F(u,v)*w(x,y) = sum(w(x,y)*[dx * d(I(x,y))/dx + dy * d(I(x,y))/dy]^2)
     Let call d(I(x,y))/dx = Ix and d(I(x,y))/dy = Iy
    (w(x,y) is window function, which can be gaussian or just a rectangular) 
    or in matrix form: F(u,v)*w(x,y)~=   [dx, dy] * M * [dx, dy]^T, where
    M = sum(w(x,y) [Ix^2 IxIy;IxIy Iy^2]) 
    From the PSA method we know, that the eigenvalues of this matrix would give us 
    the directions in which the data are mostly spread. 
    Thats why by analizing the eigenvalues of this matrix in each window, we 
    destinglish between corner, edge or none of both. 
%}
   
    %Spatial derivative calculation (Ix,Iy)
    % Create sobel operator in horizontal direction: 
    fx = [-1 0 1; -2 0 2; -1 0 1];
    % Apply it to the image
    Ix = filter2(fx,greyscale_image);
    % Create sobel operator in vertical direction:
    fy = [1 2 1; 0 0 0; -1 -2 -1];
    % Apply it to the image
    Iy = filter2(fy,greyscale_image);
    % We need to calculate also Ix^2, Iy^2, Ixy
    Ix2 = Ix.^2;
    Iy2 = Iy.^2;
    Ixy = Ix.*Iy;

    %3rd Step:  
    % as window function we choose gaussian. In the next step, we apply it to
    % the result
    h= fspecial('gaussian',[9 9], 2); 
    Ix2 = filter2(h,Ix2);
    Iy2 = filter2(h,Iy2);
    Ixy = filter2(h,Ixy);

    % Eigenvalues are very expensive to calculate, thats why we calculate a
    % parameter R(i,j) = det(M)-0.06*(trace(M))^2;
    % where 0,06 is a chosen coefficient between [0.04, 0.06]
    % Now we can destinglish between different R
    % create a matrix R of a size of the Image
    R = zeros(size(greyscale_image,1),size(greyscale_image,2));

    % set Rmax to 0 first
    Rmax = 0; 
    for i = 1:size(greyscale_image,1)
        for j = 1:size(greyscale_image,2)
            M = [Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)]; 
            R(i,j) = det(M)-0.04*(trace(M))^2;
            if R(i,j) > Rmax
                 Rmax = R(i,j);
            end
        end
    end
    % searching for local maxima as corner within the window of 3/3 
    result = zeros(size(greyscale_image,1),size(greyscale_image,2)); 
    for i = 3:size(greyscale_image,1)-2 % height 
         for j = 3:size(greyscale_image,2)-2 %width
                if R(i,j) > 0.0006*Rmax && ...
                    R(i,j)< 0.5*Rmax && ...
                    R(i,j) > R(i-1,j-1) && ...
                    R(i,j) > R(i-1,j) && ...
                    R(i,j) > R(i-1,j+1) &&...
                    R(i,j) > R(i,j-1) &&...
                    R(i,j) > R(i,j+1) && ...
                    R(i,j) > R(i+1,j-1) &&...
                    R(i,j) > R(i+1,j) && ...
                    R(i,j) > R(i+1,j+1) && ...
                    R(i,j) > R(i+2,j+2)&& ...
                    R(i,j) > R(i+2,j-2) &&...
                    R(i,j) > R(i+2,j) && ...
                    R(i,j) > R(i,j+2) && ...
                    R(i,j) > R(i,j-2) &&...
                    R(i,j) > R(i-2,j+2) &&...
                    R(i,j) > R(i-2,j-2) && ...
                    R(i,j) > R(i-2,j) 
                    result(i,j) = 1;
                end
         end
    end
    [posc, posr] = find(result == 1);
    result = [posc, posr];
end
