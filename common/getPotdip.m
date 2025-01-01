function u = getPotdip(lambda_x,lambda_y,x,y)
%GETPOTDIP(lambda_x,lambda_y,x,y) returns flow field from potential dipole


% r = bsxfun(@minus, y, x.'); 
% irr = 1./(conj(r).*r); %1/r^2
% %irr_b = 1./abs(r).^2; %same thing
% irrrr = irr.*irr;
% d1 = real(r); d2 = imag(r); 
 

c = 1/(2*pi);

N = length(lambda_x);
M = numel(y(:));
u1 = zeros(M,1);
u2 = zeros(M,1);


% for k = 1:M
%             
%     r = x-y(k); 
%     irr = 1./(conj(r).*r); %1/r^2
%     irrrr = irr.*irr;
%     d1 = real(r); 
%     d2 = imag(r); 
%     A12 = 2*d1.*d2.*irrrr.*lambda_y;
%     A21 = 2*d1.*d2.*irrrr.*lambda_x;
% 
%     u1(k) = sum((-irr+2*d1.^2.*irrrr).*lambda_x+A12);
%     u2(k) = sum((-irr+2*d2.^2.*irrrr).*lambda_y+A21);
% end

for k = 1:N
            
    r = x(k)-y; 
    irr = 1./(conj(r).*r); %1/r^2
    irrrr = irr.*irr;
    d1 = real(r); 
    d2 = imag(r); 
    A12 = 2*d1.*d2.*irrrr.*lambda_y(k);
    A21 = 2*d1.*d2.*irrrr.*lambda_x(k);

    u1 = u1 + (-irr+2*d1.^2.*irrrr).*lambda_x(k)+A12;
    u2 = u2 + (-irr+2*d2.^2.*irrrr).*lambda_y(k)+A21;
end


u = c*[u1; u2];

end