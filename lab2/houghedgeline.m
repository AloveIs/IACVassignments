function [linepar, acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose,vote_mode)
    if nargin < 8
        vote_mode = 1;
    end
    curves = extractedge(pic, scale,gradmagnthreshold,'same');
     
    magnitude = Lv(discgaussfft(pic, scale),'same');

    [linepar, acc] = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines, verbose,vote_mode);
end