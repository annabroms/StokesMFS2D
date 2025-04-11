function Lr = getLrPair(B1,B2,K)
%getLrPair help function to close the system for the unknown sources in a
% 2D Stokes mobility problem with 2-body preconditioning. 
% 
% Syntax: Lr = getLrPair(B1,B2,K) 
% 
% Description: Returns the matrix Lr for a pair of circles, constructed
% from the individual blocks Lr:=B_i*K' of each particle in the
% pair. Lr is the matrix that maps u<-lambda, i.e. it evaluates the rigid
% body velocities on the surface of the particle, u = BU, using the ansatz that U =
% -K'lambda. When we build the mobility system matrix, the sign is reversed.
%
% Note: B1,B2 and K are constructed with getKmat2D
%
% Anna Broms April 4, 2025

%There are four blocks in Lr, xx xy, yx, yy...
%Need to build projection for the pair the same way
   % Lr = blkdiag(B1*K',B2*K'); %Not the right thing! 

Lr1 = B1*K'; % xx xy; yx yy
Lr2 = B2*K';

xx1 = Lr1(1:end/2,1:end/2);
xx2 = Lr2(1:end/2,1:end/2);
yy1 = Lr1(end/2+1:end,end/2+1:end);
yy2 = Lr2(end/2+1:end,end/2+1:end);
yx1 = Lr1(end/2+1:end,1:end/2);
yx2 = Lr2(end/2+1:end,1:end/2);
xy1 = Lr1(1:end/2,end/2+1:end);
xy2 = Lr2(1:end/2,end/2+1:end);

Lrxx = [xx1 zeros(size(xx1)); zeros(size(xx1)) xx2];
Lryy = [yy1 zeros(size(xx1)); zeros(size(yy1)) yy2];
Lrxy = [xy1 zeros(size(xx1)); zeros(size(xx1)) xy2];
Lryx = [yx1 zeros(size(xx1)); zeros(size(yy1)) yx2];

Lr = [Lrxx Lrxy; Lryx Lryy];
    
end