function N1body = evaluateCoarseOnPair(q,rbase_in_c,rout_f)
%evaluateCoarseOnPair builds matrix needed in the rhs of the fine BVP to
%solve to determine pair corrections for circular particles in Stokes flow.
%Returns the map from coarse proxy Stokeslet source strengths (sorted as
%body1_x, body2_x, body1_y, body2_y) to velocities (body1_x, body2_x,
%body1_y, body2_y) evaluated in fine nodes on a particle pair. The flow
%field on one particle in the pair at the time is non-zero (the other in
%the pair) to cancel out Dirichlet data.
% 
% Syntax: N1body = evaluateCoarseOnPair(q,rbase_in_c,rout_f)
%
% Anna Broms April 4, 2025

N2 = singleLayer(rbase_in_c+q(1),rout_f(end/2+1:end),1); %source points on 1, targets on 2
N1 = singleLayer(rbase_in_c+q(2),rout_f(1:end/2),1); %source points on 2, targets on 1

%reoreder elements
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

N1body = [Nxx Nyx; Nxy Nyy];
    
end