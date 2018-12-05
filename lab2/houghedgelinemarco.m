function [linepar, acc] = ...
    houghedgelinemarco(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, voting_function)
%HOUGHEDGELINE  Extracts the edges from pic in the form of lines contained 
%               in linepar through the Hough transform

%Check validity of the input
if (nargin < 6) 
    error('Requires at least six input arguments.') 
end
if (nargin == 6)
    %Effect of the magnitude on the voting weight
    voting_function=@(x) 1;
end

%create curves
edgecurves = extractedge(pic,scale,gradmagnthreshold);
%create "valid" magnitude 
magn=Lv(discgaussfft(pic, scale),'same');
magn(1:2,:)=0;
magn(:,1:2)=0;
magn((length(pic)-1):length(pic),:)=0;
magn(:,(length(pic)-1):length(pic))=0;
%transform in line edges
[linepar, acc] = ...
    houghlinemarco(edgecurves, magn, nrho, ntheta, 0, nlines,voting_function);

end

