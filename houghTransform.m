%{
    Own implementation of a Hough transform

    Takes: A binary Image, the resolution of RHO and a vector of angles

    Returns: The accumulator image in hough space and the LUTs for both
    axes

    Author:         Anand Eichner (11808244)
%}
function [HS,THETA,RHO] = houghTransform(input,dRes,THETA)
    [sX,sY,~] = size(input);
    
    %%TODO: figure out why this implementation gives slightly different
    %%results than the builtin implementation (specifically slightly askew
    %%results)
    
    %% Calculate RHO LUT
    D = sqrt(sX*sX + sY*sY);
    q = ceil(D/dRes);
    nrho = 2*q + 1;
    RHO = linspace(-q*dRes, q*dRes, nrho);
    
    %% Calculate angles
    angles = [sind(THETA);cosd(THETA)];
    
    %% Get the positions of all true pixels of the binary image
    [pixelsX,pixelsY] = find(input);
    
    %% Calculate rho for all pixels
    distances = [pixelsX,pixelsY]*angles;
    
    %% Transform rho values to indizes
    indizes = interp1(RHO,1:nrho,distances,'nearest');
    indizes = indizes + meshgrid(0:length(THETA)-1,1:length(pixelsX)) * nrho;
    indizes = reshape(indizes,[],1);
    
    %% Accumulate all indizes
    HS = accumarray(indizes,1, [nrho*length(THETA),1],[],0);
    HS = reshape(HS, nrho, []);
    
    %figure, imshow(HS,[],'XData',THETA,'YData',RHO,'InitialMagnification','fit');
    %xlabel('\theta'), ylabel('\rho');
    %axis on, axis normal, hold on;
end