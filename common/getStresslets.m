function u = getStresslets(lambda_x,lambda_y,x,y,n1,n2)
%GETSTRESSLETS(lambda_x,lambda_y,x,y,n1,n2) returns flow field evaluated at y from
%stresslets centered at x with directions n1,n2 and strengths
%lambda_x, lambda_y. 
fmm = 0;
if fmm
    error('Not yet supported')
    [u1, u2] = stokesDLPfmm(lambda_x,lambda_y,real(x),imag(x),real(y),imag(y),n1(:),n2(:));
else
    M = numel(y(:));
    N = length(lambda_x);
    u1 = zeros(M,1);
    u2 = zeros(M,1);

% Loop over targets
%     for k = 1:M
%                         
%         rx = real(x - y(k));
%         ry = imag(x - y(k));
%         rho4 = (rx.^2 + ry.^2).^2;
%         
%         rdotq = rx.*lambda_x + ry.*lambda_y;
%         rdotn = rx.*n1 + ry.*n2;
%         
%         u1(k) = 4*sum(rdotn.*rdotq./rho4.*rx);
%         u2(k) = 4*sum(rdotn.*rdotq./rho4.*ry);
%     end

%Faster to loop over sources
    for k = 1:N
                        
        rx = real(x(k) - y);
        ry = imag(x(k) - y);
        rho4 = (rx.^2 + ry.^2).^2;
        
        rdotq = rx.*lambda_x(k) + ry.*lambda_y(k);
        rdotn = rx.*n1(k) + ry.*n2(k);
        
        u1 = u1+rdotn.*rdotq./rho4.*rx;
        u2 = u2+rdotn.*rdotq./rho4.*ry;
    end
    
    u1 = u1/pi;
    u2 = u2/pi;
end


%[u1, u2] = stokesDLPfmm(real(qwazp(:)),imag(qwazp(:)),x(:),y(:),n1(:),n2(:));
u = [-u1; -u2];

end