function Lpair = getLpair(L)
%Returns projection L for a pair of circles. 
Lxx = L(1:end/2,1:end/2);
Lxy = L(1:end/2,end/2+1:end);
Lyx = L(end/2+1:end,1:end/2);
Lyy = L(end/2+1:end,end/2+1:end);

Axx = [eye(size(Lxx))-Lxx zeros(size(Lxx)); zeros(size(Lxx)) eye(size(Lxx))-Lxx ];
Axy = [-Lxy zeros(size(Lxx)); zeros(size(Lxx)) -Lxy ];

Ayy = [eye(size(Lxx))-Lyy zeros(size(Lxx)); zeros(size(Lxx)) eye(size(Lxx))-Lyy ];
Ayx = [-Lyx zeros(size(Lxx)); zeros(size(Lxx)) -Lyx ];

Lpair = [Axx Axy; Ayx Ayy];

%debug
%     Axx2 = [Lxx zeros(size(Lxx)); zeros(size(Lxx)) Lxx ];
%     Axy2 = [Lxy zeros(size(Lxx)); zeros(size(Lxx)) Lxy ];
%     Ayy2 = [Lyy zeros(size(Lxx)); zeros(size(Lxx)) Lyy ];
%     Ayx2 = [Lyx zeros(size(Lxx)); zeros(size(Lxx)) Lyx ];
%     A2 = eye(size(Lpair))-[Axx2 Axy2; Ayx2 Ayy2]; %This is the same thing as A!

%A = [eye(size(L))-L zeros(size(L)); zeros(size(L)) eye(size(L))-L];
%%Wrong as N is ordered x,y 

end
