%{
    Start label detection.

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = label_detection(input_img, shelf_height_px, is_debug)

    img_double = im2double(input_img);
    img_grey = rgb2gray(img_double);
    
    %% APPROXIMATE LABEL DIMENSIONS / PROPERTIES
    % based on shelf height
    
    approx_label_height = shelf_height_px * 0.135;
    approx_label_width = approx_label_height * 0.583;
    min_label_height = approx_label_height * 0.65;
    max_label_height = approx_label_height * 1.2;
    min_label_width = max(20, approx_label_width * 0.2);
    max_label_width = approx_label_width * 1.2;
    
    brightness_amount_range = [7, 25];      % acceptable ratio between bright and dark areas of a label (see check_if_label())
    x_diff_range = [min_label_width, max_label_width];
    y_diff_range = [min_label_height, max_label_height];
        
    %% GLOBAL THRESHOLD
    % we apply the global threshold. all pixels lower than the th are
    % clipped to that same th value.
 
    th_img = img_grey;
    global_th = otsu_threshold(th_img);
    %th_img(th_img < global_th) = global_th;
    %global_th = otsu_threshold(th_img);
    th_mask = img_grey >= global_th;


    %% CORNER DETECTION
    % unfortunately we could not implement this on our own, because
    % the implementation of one of our colleagues was too sensitive.
    
    %img_grey_gauss = imgaussfilt(img_grey, round(val_range([3, 7], [2700, 2000], shelf_height_px)));
    %corners = corner(img_grey_gauss, 5000, ...
    %                 'FilterCoefficients', fspecial('gaussian',[round(val_range([7, 17], [2000, 2700], shelf_height_px)) 1], 2), ...
    %                 'SensitivityFactor', val_range([0.01, 0.06], [2700, 2000], shelf_height_px));
    img_grey_gaus = imgaussfilt(img_grey, 7);
    corners = corner(img_grey_gaus, 5000, 'FilterCoefficients', fspecial('gaussian',[17 1], 2), 'SensitivityFactor', 0.01);
    corners = [corners(:, 2), corners(:, 1)]; % swap columns so it's right


    %% INTEGRAL IMAGING
    % create integral image for brightness values.
    % pixels higher than the global th are assumed to be white

    img_integral_brights = generate_integral_image(double(th_mask));


    %% FIND LABELS
    
    labels = [];
    corners = sortrows(corners);            % sort corners by their y-axis

    for c = 1:size(corners, 1)
    
        neighbors = find_neighbors(corners, c, x_diff_range, y_diff_range);
        for n = 1:size(neighbors, 1)
            label_quad = create_quad(corners(c, :), neighbors(n, :));
            if check_if_label(img_grey, img_integral_brights, label_quad, brightness_amount_range)
                labels = [labels; label_quad];
            end
        end

    end
    
    %% INTERSECT & REMOVE REDUNDANT LABELS
    %{
    label_index = 1;
    while label_index < size(labels, 1)

        neighbor_index = 1;
        while neighbor_index < size(labels, 1)

            if label_index == neighbor_index
                neighbor_index = neighbor_index + 1;
                continue;
            end

            label_1 = labels(label_index, :);
            label_2 = labels(neighbor_index, :);
            if intersects_label(label_1, label_2)
                combined_label = [min([label_1(1), label_2(1)]), min([label_1(2), label_2(2)]), max([label_1(3), label_2(3)]), max([label_1(4), label_2(4)])];
                %neighbor_index = neighbor_index - 1;

                if check_for_range(combined_label(1:2), combined_label(3:4), x_diff_range, y_diff_range) && ...
                   check_if_label(img_grey, img_integral_brights, combined_label, brightness_amount_range)

                    first_index = min([label_index, neighbor_index]);
                    second_index = max([label_index, neighbor_index]);
                    labels = [labels(1:first_index-1, :); combined_label; labels(first_index+1:second_index-1, :); labels(second_index+1:end, :)];
                    neighbor_index = 1;
                    continue;
                %else
                %   labels = [labels(1:label_index-1, :); combined_label; labels(label_index+1:neighbor_index-1, :); labels(neighbor_index+1:end, :)];
                end

            end

            neighbor_index = neighbor_index + 1;
        end

        label_index = label_index + 1;
    end
    %}
    
    label_index = 1;
    while label_index < size(labels, 1)

        neighbor_index = 1;
        while neighbor_index < size(labels, 1)

            if label_index == neighbor_index
                neighbor_index = neighbor_index + 1;
                continue;
            end

            label_1 = labels(label_index, :);
            label_2 = labels(neighbor_index, :);
            if contains_label(label_1, label_2)
                combined_label = [min([label_1(1), label_2(1)]), min([label_1(2), label_2(2)]), max([label_1(3), label_2(3)]), max([label_1(4), label_2(4)])];
                
                if check_for_range(combined_label(1:2), combined_label(3:4), x_diff_range, y_diff_range) && ...
                   check_if_label(img_grey, img_integral_brights, combined_label, brightness_amount_range)

                    first_index = min([label_index, neighbor_index]);
                    second_index = max([label_index, neighbor_index]);
                    labels = [labels(1:first_index-1, :); combined_label; labels(first_index+1:second_index-1, :); labels(second_index+1:end, :)];
                    neighbor_index = 1;
                    continue;
                end

            end

            neighbor_index = neighbor_index + 1;
        end

        label_index = label_index + 1;
    end

    result = labels;
    
    %% DEBUG
    if is_debug
        corners_t = [corners(:, 2), corners(:, 1)];
        imshow(img_grey);
        hold on;
        for t = 1:size(labels, 1)
            topLeft = labels(t, 1:2);
            topRight = [labels(t, 1), labels(t, 4)];
            bottomLeft = [labels(t, 3), labels(t, 2)];
            bottomRight = labels(t, 3:4);
            plot(polyshape([topLeft(2), bottomLeft(2), bottomRight(2), topRight(2)], [topLeft(1), bottomLeft(1), bottomRight(1), topRight(1)]), 'EdgeColor', 'red', 'LineWidth', 1);
        end
        for c = 1:size(corners_t, 1)
            th = 0:pi/2:2*pi;
            r = 5;
            x = r * cos(th) + corners_t(c, 1);
            y = r * sin(th) + corners_t(c, 2);
            plot(x, y, 'r');
        end
        hold off;
    end
    
end


%% FUNCTIONS

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
function result = check_if_label(img_grey, integral_image_brights, label_quad, brights_darks_ratio_range)
    
    label_candidate = img_grey(label_quad(1):label_quad(3), label_quad(2):label_quad(4));
    dark_amount = sum(1-clear_label_border(label_candidate), 1:2);
    
    
    label_width = label_quad(4) - label_quad(2);
    label_height = label_quad(3) - label_quad(1);
    brightness_max_value  = label_width * label_height;
    
    ii_brightness_amount = integral_image_result(integral_image_brights, label_quad);    
    brightness_ratio = ii_brightness_amount / dark_amount;
    
    result = ii_brightness_amount >= brightness_max_value / 2 && ... % at least half of the label must be considered white
             brightness_ratio >= brights_darks_ratio_range(1) && ...
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
function result = clear_label_border(label_grey)

    label_binary = label_grey >= otsu_threshold(label_grey);

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
    given corners (must be top-left and bottom-right corners).

    Sources:
        http://delivery.acm.org/10.1145/810000/808600/p207-crow.pdf
        accessed on 2019/11/12

        https://en.wikipedia.org/wiki/Summed-area_table#/media/File:Summed_area_table.png
        used only for variable naming
        accessed on 2019/11/12

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = integral_image_result(integral_image, label_quad)

    corner = label_quad(1:2);
    neighbor = label_quad(3:4);

    A = integral_image(corner(1), corner(2));
    B = integral_image(corner(1), neighbor(2));
    C = integral_image(neighbor(1), corner(2));
    D = integral_image(neighbor(1), neighbor(2));
    
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
        
        MATLAB's otsuthresh-method
        accessed on 2020/01/04

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = otsu_threshold(greyscale_img)
    
    [bin_amounts, bins] = imhist(greyscale_img);
    
    bin_size = size(bins, 1);
    bin_amounts = double(bin_amounts);

    probabilities = bin_amounts / sum(bin_amounts);
    omega = cumsum(probabilities);
    mu = cumsum(probabilities .* (1:bin_size)');
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
    Checks if one quad intersects the other.
    Quads must be of style: [top-left y, top-left x, bottom-right y, bottom-right x]

    Sources:
        <ANANDS SOURCE>
        accessed on 2020/01/04

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = intersects_label(label, label_2)
    
    result = label(1) <= label_2(3) && label(2) <= label_2(4) && ...
             label(3) >= label_2(1) && label(4) >= label_2(2);

end

%{
    Checks if one quad lies completely within another.
    Works on both directions.
    Quads must be of style: [top-left y, top-left x, bottom-right y, bottom-right x]

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = contains_label(label, label_2)
    
    result = (label(1) <= label_2(1) && label(2) <= label_2(2) && ...
             label(3) >= label_2(3) && label(4) >= label_2(4)) || ...
             (label(1) >= label_2(1) && label(2) >= label_2(2) && ...
             label(3) <= label_2(3) && label(4) <= label_2(4));

end

%{
    Converts two arbitrary points to a quad, where the return value is of
    style [top-left y, top-left x, bottom-right y, bottom-right x].

    Sources:
        -

    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = create_quad(corner, neighbor)

    ii_coords_1 = [corner(1), corner(2)];
    ii_coords_12 = [corner(1), neighbor(2)];
    ii_coords_2 = [neighbor(1), neighbor(2)];
    ii_coords_21 = [neighbor(1), corner(2)];
    
    result = [min([ii_coords_1(1), ii_coords_12(1), ii_coords_2(1), ii_coords_21(1)]), ...
              min([ii_coords_1(2), ii_coords_12(2), ii_coords_2(2), ii_coords_21(2)]), ...
              max([ii_coords_1(1), ii_coords_12(1), ii_coords_2(1), ii_coords_21(1)]), ...
              max([ii_coords_1(2), ii_coords_12(2), ii_coords_2(2), ii_coords_21(2)])];
    
end