function Nout = getNPair(N1,N2)
   % N1 velocities on i, N2 maps velocities on p2

    xx1 = N1(1:end/2,1:end/2);
    xx2 = N2(1:end/2,1:end/2);
    yy1 = N1(end/2+1:end,end/2+1:end);
    yy2 = N2(end/2+1:end,end/2+1:end);
    yx1 = N1(end/2+1:end,1:end/2);
    yx2 = N2(end/2+1:end,1:end/2);
    xy1 = N1(1:end/2,end/2+1:end);
    xy2 = N2(1:end/2,end/2+1:end);
   
    Nxx = [zeros(size(xx1)) xx1;  xx2 zeros(size(xx1))];
    Nyy = [ zeros(size(xx1)) yy1; yy2 zeros(size(yy1))];
    Nxy = [ zeros(size(xx1)) xy1; xy2 zeros(size(xx1)) ];
    Nyx = [ zeros(size(xx1)) yx1; yx2 zeros(size(yy1)) ];

    Nout = [Nxx Nyx; Nxy Nyy];
    
end