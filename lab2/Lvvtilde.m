function M = Lvvtilde(image,shape)
    dx = 0.5*[-1,0,1];
    dy = 0.5*[-1;0;1];
    dxx = [1,-2,1];
    dyy = [1;-2;1];

    
    
    
    Lx = filter2(dx, image, 'same');
    Ly = filter2(dy, image, 'same');
    
    Lxy = filter2(dx, Ly, 'same');
    
    Lxx = filter2(dxx, image, 'same');
    Lyy = filter2(dyy, image, 'same');
    
    M = Lx.^2 .* Lxx + 2 * Lx .* Ly .* Lxy + Ly.^2 .* Lyy;
end