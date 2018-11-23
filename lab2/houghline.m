function [linepar acc] = houghline(curves, magnitude, nrho, ntheta, threshold, nlines, verbose)
    % Check if input appear to be valid
    
    % Allocate accumulator space
    acc = zeros(ntheta, nrho);
    img_size = size(magnitude,1);
    
    if verbose
       fprintf("Image size is %g", img_size); 
    end
    % Define a coordinate system in the accumulator space
    
    % theta -90 to 90
    % rho   - sqrt(2)* size  to sqrt(2)* size
    delta_theta = 180/ ntheta;
    delta_rho = 2 * sqrt(2)*img_size / nrho;
    rho_min = -sqrt(2)*img_size;
    
    
    % Loop over all the input curves (cf. pixelplotcurves)
    insize = size(curves, 2);
    trypointer = 1;
    numcurves = 0;
    
    while(trypointer <= insize)
        polylength = curves(2, trypointer);
        numcurves = numcurves + 1;
        trypointer = trypointer + 1;
        % For each point on each curve
        
        for polyidx = 1:polylength
            
            x = floor(curves(2, trypointer)) + 1;
            y = floor(curves(1, trypointer)) + 1;
            if verbose && x > 100
               %fprintf("\nx : %g\ty : %g",x,y); 
            end

            % Check if valid point with respect to threshold
            if(magnitude(x,y) - threshold <= 0)
                trypointer = trypointer + 1;
                continue;
            end
            
            % Optionally, keep value from magnitude image
            
            % Loop over a set of theta values
            theta_idx = 1;
            for theta = linspace(-90,90,ntheta)
                % Compute rho for each theta value
                rho = cosd(theta)*x+sind(theta)*y;
                % Compute index values in the accumulator space
                rho_idx = floor((rho-rho_min)/delta_rho)+1;
                
                
                % Update the accumulator
                acc(theta_idx, rho_idx) = acc(theta_idx, rho_idx) + 1;
                
                theta_idx = theta_idx + 1;
            end
            
            trypointer = trypointer + 1;
        end
    end
    
    if verbose
        showgrey(acc);
        return;
    end
            
    % Extract local maxima from the accumulator
    linpar = [];
    [pos, value] = locmax8(acc);
    [dummy, indexvector] = sort(value);
    nmaxima = size(value, 1);
    % Delimit the number of responses if necessary
    % Compute a line for each one of the strongest responses in the accumulator
    for idx = 1:nlines
        rhoidxacc = pos(indexvector(nmaxima - idx + 1), 1);
        thetaidxacc = pos(indexvector(nmaxima - idx + 1), 2);
        
        %compute theta
        theta = -90 + (thetaidxacc-1)* 180/ntheta;
        %compute rho
        rho = (rhoidxacc-1)*delta_rho;
        
        linpar = [linpar; rho theta];
    end
    
    
    
    % Overlay these curves on the gradient magnitude image
    % Return the output data
end

