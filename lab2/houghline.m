function [linepar acc] = houghline(curves, magnitude, nrho, ntheta, threshold, nlines, verbose,vote_mode)
    % Check if input appear to be valid
    
    % Allocate accumulator space
    acc = zeros(nrho,ntheta);
    img_size = size(magnitude,1);
    
    if verbose
       fprintf("Image size is %g", img_size); 
    end
    % Define a coordinate system in the accumulator space
    
    % theta -90 to 90
    % rho   - sqrt(2)* size  to sqrt(2)* size
    %delta_theta = 180/ ntheta;
    %delta_rho = 2 * sqrt(2)*img_size / nrho;
    rho_min = -ceil(sqrt(2)*img_size);
    
    thetas = linspace(0,180,ntheta +1);
    thetas = thetas(1:end-1);
    rhos = linspace(rho_min,-rho_min,ntheta);
    delta_rho = rhos(2) - rhos(1);
    delta_theta = thetas(2) - thetas(1);
    
    if verbose
       fprintf("drho %g \t dtheta %g", delta_rho,delta_theta); 
    end
    
    % Loop over all the input curves (cf. pixelplotcurves)
    insize = size(curves, 2);
    trypointer = 1;
    numcurves = 0;
    
    curves_matrix=pixelplotcurves(zeros(size(magnitude)),curves,1);
    
    for x=1:size(curves_matrix,2)
        for y=1:size(curves_matrix,1)
            if curves_matrix(y,x) > 0 && magnitude(y,x)>threshold 
                %update accumulator space
                
                for theta_idx = 1 : ntheta
                    theta = thetas(theta_idx);
                    % Compute rho for each theta value
                    %rho = y*cosd(theta) + x*sind(theta);
                    rho = x*cosd(theta) + y*sind(theta);
                    % Compute index values in the accumulator space
                    rho_idx = floor((rho-rho_min)/delta_rho)+1;
                    if(abs(rho-rhos(rho_idx)) >= delta_rho)
                        disp(rho-rhos(rho_idx))
                    end
                    
                    % Update the accumulator
                    switch vote_mode 

                        case 1
                            acc(rho_idx,theta_idx) = acc(rho_idx,theta_idx) + 1;
                        case 2
                            acc(rho_idx,theta_idx) = acc(rho_idx,theta_idx) + 2.7183 * log(magnitude(y,x));
                        otherwise
                            acc(rho_idx,theta_idx) = acc(rho_idx,theta_idx) + 1;
                    end
                end
            end
        end
    end
    
    
    
    if false     
        while(trypointer <= insize)
            polylength = curves(2, trypointer);
            numcurves = numcurves + 1;
            trypointer = trypointer + 1;
            % For each point on each curve

            for polyidx = 1:polylength

                x = floor(curves(2, trypointer)) + 1;
                y = floor(curves(1, trypointer)) + 1;

                % Check if valid point with respect to threshold
                if(magnitude(x,y) - threshold <= 0)
                    trypointer = trypointer + 1;
                    continue;
                end

                % Optionally, keep value from magnitude image

                % Loop over a set of theta values
                for theta_idx = 1 : ntheta
                    theta = thetas(theta_idx);
                    % Compute rho for each theta value
                    rho = cosd(theta)*x + sind(theta)*y;
                    % Compute index values in the accumulator space
                    rho_idx = floor((rho-rho_min)/delta_rho)+1;
                    if(rho-rhos(rho_idx) >= delta_rho)
                        disp(rho-rhos(rho_idx))
                    end

                    % Update the accumulator
                    acc(theta_idx, rho_idx) = acc(theta_idx, rho_idx) + 1;
                end

                trypointer = trypointer + 1;
            end
        end
    end
    
    if verbose
        showgrey(acc);
    end
    
    %%% starts here
    [pos, value] = locmax8(acc);
    [~, indexvector] = sort(value);
    nmaxima = size(value, 1);
    linepar=zeros(2,4*nlines);
    
    for i=1:nlines
        rhoidxacc = pos(indexvector(nmaxima - i + 1), 1);
        thetaidxacc = pos(indexvector(nmaxima - i + 1), 2);
        
        %compute theta
        theta = thetas(thetaidxacc); %-90 + (thetaidxacc-1)* 180/ntheta;
        %compute rho
        rho = rhos(rhoidxacc);%rho_min + (rhoidxacc-1)*delta_rho;
        
        %rho=(rho_maxima(i)-1)*deltarho-length(magnitude);
        %theta=(theta_maxima(i)-1)*deltatheta;
        x0=rho*cos(deg2rad(theta));
        y0=rho*sin(deg2rad(theta));
        dx=-sin(deg2rad(theta));
        dy=cos(deg2rad(theta));

        k = 1000;
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
    return;
    
    
    % Extract local maxima from the accumulator
    linepar = [];
    [pos, value] = locmax8(acc);
    [~, indexvector] = sort(value);
    nmaxima = size(value, 1);
    % Delimit the number of responses if necessary
    % Compute a line for each one of the strongest responses in the accumulator
    for idx = 1:nlines
        rhoidxacc = pos(indexvector(nmaxima - idx + 1), 1);
        thetaidxacc = pos(indexvector(nmaxima - idx + 1), 2);
        
        %compute theta
        theta = thetas(thetaidxacc); %-90 + (thetaidxacc-1)* 180/ntheta;
        %compute rho
        rho = rhos(rhoidxacc);%rho_min + (rhoidxacc-1)*delta_rho;
        
        linepar = [linepar; rho theta];
    end
    
    
    
    % Overlay these curves on the gradient magnitude image
    % Return the output data
end

