function [linepar, acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose)
    
    curves = extractedge(pic, scale,1,'same');
    
    disp(max(curves(1,:)));
    disp(max(curves(2,:)));
    
    deltax = 0.5*[-1,0,1];
    deltay = 0.5*[-1;0;1];
    
    dx = conv2(pic, deltax, 'same');
    dy = conv2(pic, deltay, 'same');
 
    magnitude = sqrt(dx .^2 + dy .^2);

    [linepar, acc] = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines, verbose);
end