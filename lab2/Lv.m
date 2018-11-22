function pixels = Lv(inpic,shape)
    dx = 0.5*[-1,0,1];
    dy = 0.5*[-1;0;1];

    Lx = filter2(dx, inpic, 'same');
    Ly = filter2(dy, inpic, 'same');
    pixels = Lx.^2 + Ly.^2;
end