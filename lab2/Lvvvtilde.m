function M = Lvvvtilde(image,shape)
    dx = 0.5*[-1,0,1];
    dy = 0.5*[-1;0;1];
    dxx = [1,-2,1];
    dyy = [1;-2;1];
    
    Lx = filter2(dx, image, 'same');
    Ly = filter2(dy, image, 'same');
    
    Lxx = filter2(dxx, image, 'same');
    Lyy = filter2(dyy, image, 'same');
    
    Lxxx = filter2(dx, Lxx, 'same');
    Lyyy = filter2(dy, Lyy, 'same');
    Lxxy = filter2(dxx, Ly, 'same');
    Lxyy = filter2(dx, Lyy, 'same');
    
    M = Lx.^3 .* Lxxx + 3 * Lx.^2 .* Ly .* Lxxy + 3 * Lx .* Ly.^2 .* Lxyy + Ly.^3 .* Lyyy;
end