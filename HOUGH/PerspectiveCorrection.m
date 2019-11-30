%{
    Attempts to correct the perspective of an Image based on the lines in it

    Returns: The perspective corrected image and a Struct containing an
            imRef2d object and the corners of the rectangle used to correct
            the image.

    @Author Anand Eichner
%}
function [correctedImage,spaceData] = PerspectiveCorrection(fn)
    close all;%%ffs
    
    thres = 0.4;
    peakCount = 5;
    VerticalAreaSelect = 20;


    baseImage = imread(fn);
    baseImageGray = rgb2gray(im2double(baseImage));
    
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
    edgeImageVertical = bwareafilt(edgeImageVertical,VerticalAreaSelect);
    
    %Transpose Image to reduce the houghspace to -20 to +20 degrees instead
    %of having to use a full -90 to 89 degrees
    [H,T,R] = hough(edgeImageHorizontal','RhoResolution',2,'Theta',-32.25:0.5:32);
    P  = houghpeaks(H,peakCount,'threshold',ceil(thres*max(H(:))));
    
    %Removing the zero angle seems to improve the results quite significantly
    [H2,T2,R2] = hough(edgeImageVertical,'RhoResolution',2,'Theta',-32.25:0.5:32);
    P2  = houghpeaks(H2,peakCount,'threshold',ceil(thres*max(H2(:))));
    
    %{
    figure, imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
    
    x = T(P(:,2)); y = R(P(:,1));
    plot(x,y,'s','color','white');
    %}
    
    %{
    figure, imshow(H2,[],'XData',T2,'YData',R2,'InitialMagnification','fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
    
    x2 = T2(P2(:,2)); y2 = R2(P2(:,1));
    plot(x2,y2,'s','color','white');
    %}
    
    linesHorizontal = toLines(R,T,P,true);
    linesVertical = toLines(R2,T2,P2,false);
    
    %%Plot horizontal line image
    figure, imshow(edgeImageHorizontal), hold on
    for k = 1:length(linesHorizontal)
        line = linesHorizontal(k);
        xy = findIntersections([line],[struct("P",[0,0],"D",[0,1]),struct("P",[size(edgeImageHorizontal,2),0],"D",[0,1])]);
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    end
    
    %%Plot vertical line image
    figure, imshow(edgeImageVertical), hold on
    for k = 1:length(linesVertical)
        line = linesVertical(k);
        xy = findIntersections([line],[struct("P",[0,0],"D",[1,0]),struct("P",[0,size(edgeImageHorizontal,1)],"D",[1,0])]);
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    end
    
    linesHorizontal = selectLineCandidates(linesHorizontal);
    linesVertical = selectLineCandidates(linesVertical);
    intersections = sortIntersections(findIntersections(linesVertical,linesHorizontal),true,[-1,0]);
    
    figure; imshow(baseImage), hold on;
    for i = 0:size(intersections,1)-1;
        plot(intersections([i+1,mod(i+1,4)+1],1), intersections([i+1,mod(i+1,4)+1],2),"lineWidth", 4, "color", "red");
    end
    
    %%Correct Perspective
    imageSize = size(baseImage);
    transformMatrix = fitgeotrans(...
        intersections,...
        [0 0; imageSize(2) 0; imageSize(2) imageSize(1); 0 imageSize(1)],...  %%kinda works [0 0; imageSize(1) 0; imageSize(1) imageSize(2); 0 imageSize(2)],...
        "projective");
    [correctedImage,spacialRef] = imwarp(baseImageGray,transformMatrix);
    figure(); imshow(correctedImage);
    
    spaceData.imRef = spacialRef;
    [x,y]=transformPointsForward(transformMatrix,intersections(:,1),intersections(:,2));
    spaceData.houghBounds = [x,y];
end

%{
    Takes: An array of lines

    Returns: The distances of the lines to (0,0)

    @Author Anand Eichner
%}
function [distances] = lineOriginDistances(lines)
    norms = vertcat(lines.D);
    norms = [norms(:,2),-norms(:,1)];
    points = vertcat(lines.P);
    
    distances = abs(dot(norms,points,2));
end

%{
    Takes: Hough-Transform outputs and houghpeaks

    Returns: Two best line candidates for generating a perspective plane

    @Author Anand Eichner
%}
function [lines] = selectLineCandidates(lines)
    dists = lineOriginDistances(lines);
    dists = [dists,(1:length(dists))'];
    dists = sortrows(dists,1);
    
    lines = [lines(dists(1,2)),lines(dists(length(lines),2))];
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

    Returns: A sorted list of the intersections

    @Author Anand Eichner
%}
function [intersections] = sortIntersections(intersections,reverse,zeroVector)
    if nargin > 1
        if reverse; direction = "descend"; else; direction = "ascend"; end;
    else 
        direction = "ascend";
    end
    
    if nargin > 2
        dirVec = zeroVector;
    else 
        dirVec = [0,1];
    end

    middle = sum(intersections,1)/size(intersections,1);
    dirs = intersections - middle;
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
    intersections = sortrows(intersections,3,direction);
    intersections = intersections(:,1:2);
end
