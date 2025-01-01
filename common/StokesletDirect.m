function [udirect,vdirect] = StokesletDirect(xsrc,ysrc,xtar,ytar,f1,f2,N)
%     srcEqualsTar = 0; 
     Ntarg = length(xtar); 
     udirect = zeros(Ntarg,1);
     vdirect = zeros(Ntarg,1); 
%     
%     for k = 1:Ntarg
%         
%         %if srcEqualsTar
%         %    ind = [(1:k-1) (k+1:numel(f1))];
%         %else
%           %  ind = 1:Nsrc;
%         %end
% %         rx = xtar(k) - xsrc(ind);
% %         ry = ytar(k) - ysrc(ind);
%         rx = xtar(k) - xsrc;
%         ry = ytar(k) - ysrc;
%         rho2 = rx.^2 + ry.^2;
%         rdotf = rx.*f1 + ry.*f2;
%         udirect(k) = sum(-0.5*log(rho2).*f1 + rdotf./rho2.*rx);
%         vdirect(k) = sum(-0.5*log(rho2).*f2 + rdotf./rho2.*ry);
%     end




for k = 1:N
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


