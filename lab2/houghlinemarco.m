function [linepar, acc] = ...
    houghline(curves, magnitude, nrho, ntheta, threshold, nlines, voting_function)
%HOUGHLINE      Transforms the curves in lines contained in linepar through 
%               the Hough transform

%Check validity of the input
if (nargin < 6) 
    error('Requires at least six input arguments.') 
end
if (nargin == 6)
    %Effect of the magnitude on the voting weight
    voting_function=@(x) 1;
end

%Allocate accumulator space
acc=zeros(nrho,ntheta);

%Define a coordinate system in the accumulator space
%In the worst case scenario rho goes from -x(y)_max to sqrt(2)*x(y)_max
%I want a fixed deltarho, so I divide all the dynamical range in advance
%(it's always possible to call the function with a smaller/bigger nrho)
deltarho=length(magnitude)*(sqrt(2)+1)/nrho;
thetamax=180;   %deg
deltatheta=thetamax/ntheta; %theta goes from 0 to (ntheta-1)*thetamax


%Transform curve into a mask
curves_matrix=pixelplotcurves(zeros(size(magnitude)),curves,1);

%Summarize all the information in one matrix
curves_weighted=curves_matrix.*magnitude;

%Loop over all input curves
for i=1:size(curves_matrix,1)
    for j=1:size(curves_matrix,2)
        if curves_weighted(i,j)>threshold
            %update accumulator space
            vote_weight= 1;
            for h=1:ntheta
                rho=j*cos(deg2rad((h-1)*deltatheta)) + i*sin(deg2rad((h-1)*deltatheta));
                index_rho=ceil((rho+length(magnitude))/deltarho);
                acc(index_rho,h)=acc(index_rho,h) + vote_weight;
            end
        end
    end
end

%Extract local maxima
%acc=binsepsmoothiter(acc, 0.2, 1);
rho_maxima=zeros(nlines,1);
theta_maxima=zeros(nlines,1);
[pos, value] = locmax8(acc);
[dummy, indexvector] = sort(value);
nmaxima = size(value,1);

for i=1:nlines
    rho_maxima(i) = pos(indexvector(nmaxima - i +1), 1);
    theta_maxima(i) = pos(indexvector(nmaxima - i +1), 2);
end

%Create lines in the curves format
linepar=zeros(2,4*nlines);
for i=1:nlines
    rho=(rho_maxima(i)-1)*deltarho-length(magnitude);
    theta=(theta_maxima(i)-1)*deltatheta;
    x0=rho*cos(deg2rad(theta));
    y0=rho*sin(deg2rad(theta));
    dx=-sin(deg2rad(theta));
    dy=cos(deg2rad(theta));
    
    if dx*dy~=0
        k=ceil(max(abs(length(magnitude)/dx),abs(length(magnitude)/dy)));
    else 
        k=length(magnitude);
    end
    dx=dx*k;
    dy=dy*k;
    
    linepar(1, 4*(i-1) + 1) = 0; % level, not significant
    linepar(2, 4*(i-1) + 1) = 3; % number of points in the curve
    linepar(1, 4*(i-1) + 2) = x0-dx;
    linepar(2, 4*(i-1) + 2) = y0-dy;
    linepar(1, 4*(i-1) + 3) = x0;
    linepar(2, 4*(i-1) + 3) = y0;
    linepar(1, 4*(i-1) + 4) = x0+dx;
    linepar(2, 4*(i-1) + 4) = y0+dy;
end

end