%{
    Attempts to correct the perspective of an Image based on the lines in it

    Returns: The perspective corrected image and a Struct containing an
            imRef2d object, the corners of the rectangle used to correct
            the image, and all horizontal Lines that were used(found).

    @Author Anand Eichner
%}
function [correctedImage,spaceData] = PerspectiveCorrection(baseImage)
    baseImageGray = rgb2gray(im2double(baseImage));
    
    thres = 0.4;
    peakCount = 5;
    VerticalAreaSelect = 20;
    
    %%ONLY change this if you know what it is for (this is a safety feature)
    %%should be around 1.5 to 2 in most cases, 5 in extreme cases
    IMAGE_MAX_REL_SIZE = 3;
    %%reduces overall output size
    ENABLE_PRE_CORR_UPSCALING = true;
    %%Set to true to display different stages of the process on screen
    DEBUG = true;
    if(DEBUG); close all; end%ffs
    
    %% Preprocessing
    
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
    
    %% Hough Transform
    %Transpose Image to reduce the houghspace to -20 to +20 degrees instead
    %of having to use a full -90 to 89 degrees
    [H,T,R] = hough(edgeImageHorizontal','RhoResolution',2,'Theta',-32.25:0.5:32);
    P  = houghpeaks(H,peakCount,'threshold',ceil(thres*max(H(:))));
    
    %Removing the zero angle seems to improve the results quite significantly
    [H2,T2,R2] = hough(edgeImageVertical,'RhoResolution',2,'Theta',-32.25:0.5:32);
    P2  = houghpeaks(H2,peakCount,'threshold',ceil(thres*max(H2(:))));
    
    %% Find a rectangle in the plane of the shelf
    linesHorizontal = toLines(R,T,P,true);
    linesVertical = toLines(R2,T2,P2,false);
    verticalBounds = selectLineCandidates(linesHorizontal);
    horizontalBounds = selectLineCandidates(linesVertical);
    rectangleBounds = sortIntersections(findIntersections(horizontalBounds,verticalBounds),true,[-1,0]);
    
    %%Optionally upscale bounds to fit the image as good as possible
    imageSize = size(baseImage);
    if (ENABLE_PRE_CORR_UPSCALING); scaledBounds = upscalePerspectiveRectangle(rectangleBounds,imageSize);
    else; scaledBounds = rectangleBounds;
    end
    
    %% Correct Perspective
    transformMatrix = fitgeotrans(...
        scaledBounds,...
        [0 0; imageSize(2) 0; imageSize(2) imageSize(1); 0 imageSize(1)],...  %%kinda works [0 0; imageSize(1) 0; imageSize(1) imageSize(2); 0 imageSize(2)],...
        "projective");
    %% TODO: May need to do some aspect ratio correction
    
    %%Prevent Matlab from locking up your PC if it finds an invalid
    %%rectangle (too small)
    [xlim, ylim] = outputLimits(transformMatrix,[1 imageSize(2)],[1 imageSize(1)]);
    if(((xlim(2)-xlim(1)) * (ylim(2)-ylim(1))) / (imageSize(2)*imageSize(1)) > IMAGE_MAX_REL_SIZE*IMAGE_MAX_REL_SIZE)
        correctedImage = baseImage;
        spacialRef = imref2d();
    else
        %% Warp image using the transformmatrix
        [correctedImage,spacialRef] = imwarp(baseImage,transformMatrix);
    end
    
    %% Write additional return data
    spaceData.imRef = spacialRef;
    imOffset = [-spacialRef.XWorldLimits(1),-spacialRef.YWorldLimits(1)];
    spaceData.originalBounds = [0 0; imageSize(2) 0; imageSize(2) imageSize(1); 0 imageSize(1)] + repmat(imOffset,4,1);
    spaceData.horizontals = transformLines(linesHorizontal,transformMatrix,imOffset);
    
    %% Debug visualizations (These are completely post work)
    if DEBUG
        figure, imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
        xlabel('\theta'), ylabel('\rho');
        axis on, axis normal, hold on;
        x = T(P(:,2)); y = R(P(:,1));
        plot(x,y,'s','color','white');

        figure, imshow(H2,[],'XData',T2,'YData',R2,'InitialMagnification','fit');
        xlabel('\theta'), ylabel('\rho');
        axis on, axis normal, hold on;
        x2 = T2(P2(:,2)); y2 = R2(P2(:,1));
        plot(x2,y2,'s','color','white');

        figure, imshow(edgeImageHorizontal), hold on
        for k = 1:length(linesHorizontal)
            line = linesHorizontal(k);
            xy = findIntersections([line],[struct("P",[0,0],"D",[0,1]),struct("P",[size(edgeImageHorizontal,2),0],"D",[0,1])]);
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        end

        figure, imshow(edgeImageVertical), hold on
        for k = 1:length(linesVertical)
            line = linesVertical(k);
            xy = findIntersections([line],[struct("P",[0,0],"D",[1,0]),struct("P",[0,size(edgeImageHorizontal,1)],"D",[1,0])]);
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        end

        figure; imshow(baseImage), hold on;
        for i = 0:size(rectangleBounds,1)-1
            plot(rectangleBounds([i+1,mod(i+1,4)+1],1), rectangleBounds([i+1,mod(i+1,4)+1],2),"lineWidth", 4, "color", "red");
        end
        
        figure(); imshow(correctedImage), hold on;
        for k = 1:length(spaceData.horizontals)
            line = spaceData.horizontals(k);
            xy = findIntersections([line],[struct("P",[0,0],"D",[0,1]),struct("P",[size(correctedImage,2),0],"D",[0,1])]);
            plot(xy(:,1),xy(:,2),'LineWidth',3,'Color','green');
        end
        for i = 0:size(spaceData.originalBounds,1)-1
            plot(spaceData.originalBounds([i+1,mod(i+1,4)+1],1), spaceData.originalBounds([i+1,mod(i+1,4)+1],2),"lineWidth", 2, "color", "red");
        end
    end
end

%% Local Function Definitions

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
    Takes: An array of lines

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
    if transposed; a = 90-a; end
    nd = R(P(:,1));
    if transposed; nd = -nd; end
    d = [sind(a'),cosd(a')];
    p = repmat(nd',1,2) .* [d(:,2),-d(:,1)];
    lines = [];
    for i = 1:size(p,1)
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
function [intersection,t] = findIntersection(lineA,lineB)
    d=lineA.D(1)*lineB.D(2)-lineA.D(2)*lineB.D(1);
    if d == 0; intersection = []; t=[]; return; end
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
    %% Process arguments
    if nargin > 1
        if reverse; direction = "descend"; else; direction = "ascend"; end
    else 
        direction = "ascend";
    end
    
    if nargin > 2
        dirVec = zeroVector;
    else 
        dirVec = [0,1];
    end

    %% Calculate weight point and the direction vectors of all points relative to it
    middle = sum(intersections,1)/size(intersections,1);
    dirs = intersections - middle;
    dirssqrd = dirs.*dirs;
    dirs = dirs ./ repmat(sqrt(dirssqrd(:,1)+dirssqrd(:,2)),1,2);
    
    %% Compute angles
    c = cross([dirs,zeros(length(dirs),1)],repmat([dirVec,0],length(dirs),1),2);
    c = c(:,3);
    d = dot(dirs,repmat(dirVec,length(dirs),1),2);
    angle = acosd(d);
    angle = abs(360*(c<0) - angle);
    
    %% Sort by angles
    intersections = [intersections,angle];
    intersections = sortrows(intersections,3,direction);
    intersections = intersections(:,1:2);
end

%{
    Takes: An array of lines and a transformmatrix

    Returns: The transformed lines

    @Author Anand
%}
function [lines] = transformLines(lines,tform,imOffset)
    points = vertcat(lines.P);
    points2 = vertcat(lines.D) + points;
    
    points = transformPointsForward(tform,points);
    points2 = transformPointsForward(tform,points2);
    
    dirs = points2 - points;
    dirs = dirs ./ repmat(sqrt(dot(dirs,dirs,2)),1,2);
    points = points + repmat(imOffset,length(lines),1);
    lines = [];
    for i = 1:size(points,1)
        lines = [lines,struct("P",points(i,:),"D",dirs(i,:))];
    end
end

%{
    Takes: A line and a box (array of four lines)

    Returns: The two points of intersection

    @Author Anand Eichner
%}
function [intersections] = intersectBox(line,box)
    intersections = [];
    for j = 1:size(box,2)
        [intersection,t] = findIntersection(box(j),line);
        intersections = [intersections; intersection, t];
    end
    intersections = sortrows(intersections,3);
    intersections = intersections(2:3,1:2);
end

%{
    Takes: Two Points

    Returns: A line going through both points in parameter form

    @Author Anand
%}
function [lines] = pointsToLine(pointA,pointB)
    dirs = pointB - pointA;
    dirs = dirs ./ repmat(sqrt(dot(dirs,dirs,2)),1,2);
    lines = [];
    for i = 1:size(pointA,1)
        lines = [lines,struct("P",pointA(i,:),"D",dirs(i,:))];
    end
end

%{
    Takes: A point and the size vector of the image

    Returns: The sector that the point is in

    Notes: 0  1  2
           7 -1  3
           6  5  4

    @Author Anand
%}
function [sector] = pointSector(point,imageSize)
    if(point(1)<0)
        if(point(2)<0)
            sector = 0;
        elseif(point(2)>imageSize(1))
            sector = 2;
        else
           sector = 1; 
        end
    elseif(point(1)>imageSize(2))
        if(point(2)<0)
            sector = 6;
        elseif(point(2)>imageSize(1))
            sector = 4;
        else
           sector = 5; 
        end
    else
        if(point(2)<0)
            sector = 7;
        elseif(point(2)>imageSize(1))
            sector = 3;
        else
           sector = -1; 
        end
    end
end

%{
    Takes: matrix with two points and a point

    Returns: the other point

    @Author Anand
%}
function [ret] = selectOtherPoint(points,point)
    if(points(1,:) == point)
        ret = points(2,:);
    else
        ret = points(1,:);
    end
end

%{
    Test Data for the upscale Function
    %upscalePerspectiveRectangle([0.22,0.3;0.78,0.42;0.86,0.68;0.27,0.6],[1,1]) %% case 4
    %upscalePerspectiveRectangle([0.22,0.3;0.78,0.42;0.83,0.63;0.27,0.6],[1,1]) %% case 2
    %upscalePerspectiveRectangle([0.25,0.34;0.76,0.33;0.69,0.55;0.22,0.53],[1,1]) %% case 2.2
    %upscalePerspectiveRectangle([0.36,0.39;0.58,0.43;0.56,0.59;0.34,0.55],[1,1]) %% case 3
    %upscalePerspectiveRectangle([0.36,0.39;0.61,0.35;0.67,0.53;0.42,0.58],[1,1]) %% case 3.2
    %upscalePerspectiveRectangle([0.35,0.32;0.74,0.3;0.78,0.57;0.32,0.65],[1,1]) %% case 1
%}
%{
    Takes: An a sorted list of four intersections defining the bounds of a
    rectange in space and the size of the image

    Returns: The the bounds scaled up to fit the ImageSize while keeping
    the same perspective distortion

    @Author Anand
%}
function [outbounds] = upscalePerspectiveRectangle(bounds,imageSize)

    %%Set to true to display different stages of the process on screen
    DEBUG = true;

    imageBoundingBox = [...
        struct("P",[0,0],"D",[0,1]),...
        struct("P",[imageSize(2),0],"D",[0,1]),...
        struct("P",[0,0],"D",[1,0]),...
        struct("P",[0,imageSize(1)],"D",[1,0])
        ];
    
    imageRect = [0,0;imageSize(2),0;imageSize(2),imageSize(1);0,imageSize(1)];
    
    %%Decide on a new strategy for bounds that already extend beyond the imageSize
    %%for now just skip the scaling
    for i = 1:size(bounds,1)
        point = bounds(i,:);
        if(point(1)<0 || point(1)>imageSize(2) || point(2)<0 || point(2)>imageSize(1))
            outbounds = bounds;
            return;
        end
    end
    
    %% create lines going along all directions of the quad
    topLine = pointsToLine(bounds(1,:),bounds(2,:));
    rightLine = pointsToLine(bounds(2,:),bounds(3,:));
    bottomLine = pointsToLine(bounds(3,:),bounds(4,:));
    leftLine = pointsToLine(bounds(4,:),bounds(1,:));
    
    %% find escape Points of the rectangle
    horizontalEscape = findIntersection(topLine,bottomLine);
    verticalEscape = findIntersection(leftLine,rightLine);
    
    heSector = pointSector(horizontalEscape,imageSize);
    veSector = pointSector(verticalEscape,imageSize);
    
    %% select Strategy
    if(heSector==-1 || veSector==-1) %% This is a very rare and very complex case that I am not going to do right now
        outbounds = bounds;
    elseif(mod(heSector,2)==1 && mod(veSector,2)==1)
        if(mod(heSector,4) ~= mod(veSector,4)) %% One internal point symmetrical approach
            pointA = imageRect(mod((heSector+veSector)/4,4)+1,:);
            pointB = selectOtherPoint(intersectBox(pointsToLine(pointA,verticalEscape),imageBoundingBox),pointA);
            pointC = selectOtherPoint(intersectBox(pointsToLine(pointA,horizontalEscape),imageBoundingBox),pointA);
            pointD = findIntersection(pointsToLine(pointB,horizontalEscape),pointsToLine(pointC,verticalEscape));
            outbounds = sortIntersections([pointA;pointB;pointC;pointD],true,[-1,0]);
        else %% Edge case
            outbounds = bounds;
        end
    elseif(mod(heSector,2)==1 && mod(veSector,2)==0 ||... %% One internal point asymmetrical approach
            mod(heSector,2)==0 && mod(veSector,2)==1)
        %% classify corner and edge sector
        if (mod(heSector,2)==0)
            corner = heSector;
            edge = veSector;
        else 
            corner = veSector;
            edge = heSector;
        end
            
        %% Select correct starting point
        if(mod(abs(heSector-veSector),6)==1)
            pointA = imageRect(max(mod(fix(heSector/2),4),mod(fix(veSector/2),4))+1,:);
        else
            pointA = imageRect(mod(corner/2,4)+1,:);
        end
        direction = sign(corner-edge);
        
        %% Select order to process thing in
        if (xor(corner==2 || corner == 6,direction==-1))
            firstEscape = verticalEscape;
            secondEscape = horizontalEscape;
        else 
            firstEscape = horizontalEscape;
            secondEscape = verticalEscape;
        end
        
        pointB = selectOtherPoint(intersectBox(pointsToLine(pointA,firstEscape),imageBoundingBox),pointA);
        pointC = selectOtherPoint(intersectBox(pointsToLine(pointB,secondEscape),imageBoundingBox),pointB);
        pointD = findIntersection(pointsToLine(pointA,secondEscape),pointsToLine(pointC,firstEscape));
        
        outbounds = sortIntersections([pointA;pointB;pointC;pointD],true,[-1,0]);
        
    elseif(mod(heSector,2)==0 && mod(veSector,2)==0)
        if(mod(heSector,4) == mod(veSector,4)) %% Two internal points diagonal approach
            pointA = imageRect(heSector/2+1,:);
            pointB = imageRect(mod(heSector/2+2,4)+1,:);
            pointC = findIntersection(pointsToLine(pointA,horizontalEscape),pointsToLine(pointB,verticalEscape));
            pointD = findIntersection(pointsToLine(pointB,horizontalEscape),pointsToLine(pointA,verticalEscape));
            outbounds = sortIntersections([pointA;pointB;pointC;pointD],true,[-1,0]);
        else %% No internal points general approach
            lineA = topLine;
            pointA = intersectBox(lineA,imageBoundingBox);
            pointA = pointA(2,:);
            
            pointB = selectOtherPoint(intersectBox(pointsToLine(pointA,verticalEscape),imageBoundingBox),pointA);
            pointC = selectOtherPoint(intersectBox(pointsToLine(pointB,horizontalEscape),imageBoundingBox),pointB);
            pointD = findIntersection(lineA,pointsToLine(pointC,verticalEscape));

            outbounds = sortIntersections([pointA;pointB;pointC;pointD],true,[-1,0]);
        end
    else
        outbounds = bounds; %% Also a very rare and very complex case that I am not going to do right now
    end
    
    %% plot stuff for debugging
    if(DEBUG)
        figure, set(gca, 'YDir','reverse'), hold on;
        plot([0,imageSize(2)],[0,0],'LineWidth',3,'Color','green');
        plot([0,imageSize(2)],[imageSize(1),imageSize(1)],'LineWidth',3,'Color','green');
        plot([imageSize(2),imageSize(2)],[imageSize(1),0],'LineWidth',3,'Color','green');
        plot([0,0],[imageSize(1),0],'LineWidth',3,'Color','green');

        plot(bounds(1:2,1),bounds(1:2,2),'LineWidth',3,'Color','red');
        plot(bounds(2:3,1),bounds(2:3,2),'LineWidth',3,'Color','red');
        plot(bounds(3:4,1),bounds(3:4,2),'LineWidth',3,'Color','red');
        plot(bounds([4,1],1),bounds([4,1],2),'LineWidth',3,'Color','red');

        xy = intersectBox(topLine,imageBoundingBox);
        xy2 = intersectBox(leftLine,imageBoundingBox);
        
        linesA = pointsToLine(xy,repmat(verticalEscape,2,1));
        linesB = pointsToLine(xy2,repmat(horizontalEscape,2,1));
        outerBounds = sortIntersections(findIntersections(linesA,linesB),true,[-1,0]);
        
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
        plot(xy2(:,1),xy2(:,2),'LineWidth',2,'Color','blue');
        for i = 0:size(outerBounds,1)-1
            plot(outerBounds([i+1,mod(i+1,4)+1],1), outerBounds([i+1,mod(i+1,4)+1],2),"lineWidth", 1, "color", "red");
        end
        
        for i = 0:size(outbounds,1)-1
            plot(outbounds([i+1,mod(i+1,4)+1],1), outbounds([i+1,mod(i+1,4)+1],2),"lineWidth", 1, "color", "green");
        end
    end
end