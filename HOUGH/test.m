%{
    Attempts to correct the perspective of an Image based on the lines in it
    @Author Anand Eichner
%}
function [] = testHough(fn,thres,peakCount)
    close all;%%ffs

    baseImage = imread(fn);
    baseImageGray = rgb2gray(im2double(baseImage));
    
    linesA = [struct("P",[0,0],"D",[1,0]),struct("P",[0,1],"D",[1,0])];
    linesB = [struct("P",[0,0],"D",[0,1]),struct("P",[1,0],"D",[0,1])];
    a = findIntersections(linesA,linesB)
    sortIntersections(a)
    
    %Only take horizontal Edges
    %make sure to extend edges and not zero pad them to avoid detecting
    %edges along the image borders
    edgeImageHorizontal = imbinarize(...
        abs(...
        imfilter(baseImageGray,[1,2,1;0,0,0;-1,-2,-1],"replicate")...
        )...
        );
    edgeImageHorizontal = bwareaopen(edgeImageHorizontal,100);
    
    %Only take vertical Edges
    %make sure to extend edges and not zero pad them to avoid detecting
    %edges along the image borders
    edgeImageVertical = imbinarize(...
        abs(...
        imfilter(baseImageGray,[1,2,1;0,0,0;-1,-2,-1]',"replicate")...
        )...
        );
    edgeImageVertical = bwareaopen(edgeImageVertical,500);
    edgeImageVertical = bwareafilt(edgeImageVertical,peakCount);
    
    %Transpose Image to reduce the houghspace to -20 to +20 degrees instead
    %of having to use a full -90 to 89 degrees
    [H,T,R] = hough(edgeImageHorizontal','RhoResolution',2,'Theta',-32.25:0.5:32);
    
    %Removing the zero angle seems to improve the results quite significantly
    [H2,T2,R2] = hough(edgeImageVertical,'RhoResolution',2,'Theta',-32.25:0.5:32);
    
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
    
    %%Plot horizontal line image
    figure, imshow(edgeImageHorizontal), hold on
    linesHorizontal = toLines(R,T,P,true);
    for k = 1:length(linesHorizontal)
        line = linesHorizontal(k);
        xy = findIntersections([line],[struct("P",[0,0],"D",[0,1]),struct("P",[size(edgeImageHorizontal,2),0],"D",[0,1])]);
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    end
    
    %%Plot vertical line image
    figure, imshow(edgeImageVertical), hold on
    linesVertical = toLines(R2,T2,P2,false);
    for k = 1:length(linesVertical)
        line = linesVertical(k);
        xy = findIntersections([line],[struct("P",[0,0],"D",[1,0]),struct("P",[0,size(edgeImageHorizontal,1)],"D",[1,0])]);
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    end
    
    %%Correct Perspective
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

%{
    Takes: Hough-Transform outputs and houghpeaks

    Returns: An array of Lines

    @Author Anand Eichner
%}
function [lines] = toLines(R,T,P,transposed)
    a = -T(P(:,2));
    if transposed a = 90-a; end
    nd = R(P(:,1));
    if transposed nd = -nd; end
    d = [sind(a'),cosd(a')];
    p = repmat(nd',1,2) .* [d(:,2),-d(:,1)];
    lines = [];
    for i = 1:size(p,1);
        lines = [lines,struct("P",p(i,:),"D",d(i,:))];
    end
end

%{
    Takes: Two arrays of lines

    Returns: An array of all Intersections between the two groups

    @Author Anand Eichner
%}
function [intersections] = findIntersections(lineGroupA,lineGroupB)
    intersections = [];
    for i = 1:size(lineGroupA,2)
        for j = 1:size(lineGroupB,2)
            intersections = [intersections; findIntersection(lineGroupA(i),lineGroupB(j))];
        end
    end
end

%{
    Takes: Two Lines

    Returns: The intersection Point
        (an empty Vector if coliniear or parallel)

    @Author Anand Eichner
%}
function [intersection] = findIntersection(lineA,lineB)
    d=lineA.D(1)*lineB.D(2)-lineA.D(2)*lineB.D(1);
    if d == 0 intersection = []; return; end
    t = (lineB.P(1) - lineA.P(1))*lineB.D(2) - (lineB.P(2) - lineA.P(2))*lineB.D(1);
    t = t/d;
    intersection = lineA.P + t * lineA.D;
end

%{
    Takes: A list of intersections
        (optionally a directionVector at which to start)

    Returns: A sorted list of the intersections (clockwise)

    @Author Anand Eichner
%}
function [intersections] = sortIntersections(intersections,zeroVector)
    if nargin > 2
      dirVec = zeroVector;
    else
      dirVec = [0,1];
    end

    middle = sum(intersections,1)/size(intersections,1);
    dirs = intersections - middle;
    %%normalize directions
    dirssqrd = dirs.*dirs;
    dirs = dirs ./ repmat(sqrt(dirssqrd(:,1)+dirssqrd(:,2)),1,2);
    
    %%Compute angles
    c = cross([dirs,zeros(length(dirs),1)],repmat([dirVec,0],length(dirs),1),2);
    c = c(:,3);
    d = dot(dirs,repmat(dirVec,length(dirs),1),2);
    angle = acosd(d);
    angle = abs(360*(c<0) - angle);
    
    %%Sort by angles
    intersections = [intersections,angle];
    intersections = sortrows(intersections,3);
    intersections = intersections(:,1:2)
end