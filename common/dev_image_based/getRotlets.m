function [uR] = getRotlets(tau,x,y)
%GETROTLETS(tau,x,y) computes 2D rotlet flow field in target points y due to
%scalar torques tau at source locations x

uR = zeros(length(y),1);

% Loop over sources
for i = 1:length(tau)
    r = y - x(i);

    %Anna: 2023-01-09: %To follow the convention for the perp operator: (-f2,f1)
    uR = uR - (imag(r) - 1i*real(r))*tau(i)./abs(r).^2/(4*pi);
        
end

uR = [real(uR); imag(uR)];


end

