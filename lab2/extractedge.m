function edgecurves = extractedge(inpic, scale, threshold, shape)
    thresh = true;
    if(nargin < 4)
        thresh = false;
        shape = 'same';
    end
    
    inpic = discgaussfft(inpic, scale);
    
    curves = zerocrosscurves(Lvvtilde(inpic,shape), -Lvvvtilde(inpic,shape));
    
    if(thresh)
        curves = thresholdcurves(curves, (Lv(inpic,shape)-threshold));
    end
    
    edgecurves = curves;
end