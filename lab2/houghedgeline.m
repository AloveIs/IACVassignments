function [linepar acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose)
    
    curves = extractedge(pic, scale);
    magnitude;

    [linepar acc] = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines, verbose);
end