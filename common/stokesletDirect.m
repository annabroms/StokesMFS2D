function [udirect,vdirect] = stokesletDirect(xsrc,ysrc,xtar,ytar,f1,f2,N)
%STOKESLETDIRECT(xsrc,ysrc,xtar,ytar,f1,f2,N) computes contribution from
%N sources located at xsrc, ysrc, with strengths fx, fy in specified target points xtar,ytar.
%Returns flow udirect,vdirect, the x and y components of the flow. Sources
%and targets assumed to be different points.

%Todo: write to support multi d source

assert(numel(xsrc)==N,'Wrong number of x sources');
assert(numel(ysrc)==N,'Wrong number of y sources');
assert(numel(f1)==N,'Wrong number of x source strengths');
assert(numel(f2)==N,'Wrong number of y source strengths');

udirect = zeros(size(ytar));
vdirect = zeros(size(ytar)); 


for k = 1:N % loop over sources
     rx = xtar - xsrc(k);
     ry = ytar - ysrc(k);
     rho2 = rx.^2 + ry.^2;
     rdotf = rx.*f1(k) + ry.*f2(k);
     logrho2 = log(rho2);
     udirect = udirect-0.5*logrho2.*f1(k) + rdotf./rho2.*rx;
     vdirect = vdirect-0.5*logrho2.*f2(k) + rdotf./rho2.*ry;
end


    vdirect = vdirect/4/pi;
    udirect = udirect/4/pi; 
end

