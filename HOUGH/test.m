%{
    Test for the hough transforms, edge detection and perspective
    correction
%}

function [] = testHough(fn,thres,peakCount)
    baseImage = imread(fn);
    baseImageGray = rgb2gray(im2double(baseImage));
    %edgeImage = edge(baseImageGray,"canny",[0.1,0.38]);
    
    %Only take horizontal Edges
    %make sure to extend edges and not zero pad them to avoid detecting
    %edges along the image borders
    edgeImageHorizontal = bwareaopen(imbinarize(abs(imfilter(baseImageGray,[1,2,1;0,0,0;-1,-2,-1],"replicate"))),100);
    edgeImageVertical = bwareaopen(imbinarize(abs(imfilter(baseImageGray,[1,2,1;0,0,0;-1,-2,-1]',"replicate"))),500);

    %Transpose Image to reduce the houghspace to -20 to +20 degrees instead
    %of having to use a full -90 to 89 degrees
    [H,T,R] = hough(edgeImageHorizontal','RhoResolution',0.5,'Theta',-32:0.5:32);
    
    [H2,T2,R2] = hough(edgeImageHorizontal,'RhoResolution',0.5,'Theta',-32:0.5:32);
    
    figure, imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;

    P  = houghpeaks(H,peakCount,'threshold',ceil(thres*max(H(:))));
    x = T(P(:,2)); y = R(P(:,1));
    plot(x,y,'s','color','white');
    
    figure, imshow(H2,[],'XData',T2,'YData',R2,'InitialMagnification','fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
    
    P2  = houghpeaks(H2,peakCount,'threshold',ceil(thres*max(H2(:))));
    x2 = T2(P2(:,2)); y2 = R2(P2(:,1));
    plot(x2,y2,'s','color','white');

    lines2 = houghlines(edgeImageVertical,T2,R2,P2,'FillGap',50,'MinLength',150);
    
    lines = houghlines(edgeImageHorizontal',T,R,P,'FillGap',50,'MinLength',150);
    
    figure, imshow(edgeImageVertical), hold on
    max_len = 0;
    for k = 1:length(lines2)
       xy = [lines2(k).point1; lines2(k).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines2(k).point1 - lines2(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    
    figure, imshow(edgeImageHorizontal), hold on
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,2),xy(:,1),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(2,2),xy(2,1),'x','LineWidth',2,'Color','yellow');
       plot(xy(1,2),xy(1,1),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    
    transformMatrix = fitgeotrans(...
        [0 0; 0 1; 1 1; 1 0],...
        [0 0; 0 1; 0.8 4; 0.8 -1],...
        "projective");
    tranformCoordsys = imref2d(size(baseImage),[0,1],[0,1]);
    tranformImageOutputBounds=imref2d(size(baseImageGray),[1 size(baseImageGray,2)],[1 size(baseImageGray,1)]);
    transformedImage = imwarp(baseImageGray,tranformCoordsys,transformMatrix); %"OutputView", tranformImageOutputBounds
    figure();
    imshow(transformedImage);
end