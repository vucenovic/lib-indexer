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
