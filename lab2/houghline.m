function [linepar acc] = houghline(curves, magnitude, nrho, ntheta, threshold, nlines, verbose)
    % Check if input appear to be valid
    
    % Allocate accumulator space
    acc = zeros(ntheta, nrho);
    % Define a coordinate system in the accumulator space
    
    % theta -90 to 90
    % rho   - sqrt(2)* size  to sqrt(2)* size
    
    delta_rho = 2*sqrt(2)*img_size / nrho;
    
    
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
            
            
            x = curves(2, trypointer);
            y = curves(1, trypointer);
            
            % Check if valid point with respect to threshold
            if(magnitude(x,y) - threshold <= 0)
                continue;
            end
            
            % Optionally, keep value from magnitude image
            
            % Loop over a set of theta values
            theta_idx = 1;
            for theta = linspace(-90,90,ntheta)
                % Compute rho for each theta value
                rho = cos(theta)*x+sin(theta)*y;
                % Compute index values in the accumulator space
                rho_idx = round(rho/delta_rho,0);
                
                % Update the accumulator
                acc(theta_idx, rho_idx) = acc(theta_idx, rho_idx) + 1;
                
                theta_idx = theta_idx + 1;
            end
            
            trypointer = trypointer + 1;
        end
    end
        
            
    % Extract local maxima from the accumulator
    
    [pos value] = locmax8(acc);
    [dummy indexvector] = sort(value);
    nmaxima = size(value, 1);
    % Delimit the number of responses if necessary
    % Compute a line for each one of the strongest responses in the accumulator
    for idx = 1:nlines
        rhoidxacc = pos(indexvector(nmaxima - idx + 1), 1);
        thetaidxacc = pos(indexvector(nmaxima - idx + 1), 2);
        
        %compute theta
        theta = -90 + (rhoidxacc-1)* 180/nrho;
        %compute rho
        
        
    end
    
    
    
    % Overlay these curves on the gradient magnitude image
    % Return the output data
end

