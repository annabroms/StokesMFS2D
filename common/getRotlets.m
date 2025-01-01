function [uR] = getRotlets(lambda,x,y)
%Rotlet contribution

uR = zeros(length(y),1);

for i = 1:length(zsrc)
    r = y - x(i);

    %Anna: 2023-01-09: %To follow the convention for the perp operator: (-f2,f1)
    uR = uR - (imag(r) - 1i*real(r))*lambda(i)./abs(r).^2/(4*pi);
    
    %David: uR = uR + (-imag(r) + 1i*real(r))*torques(i)./abs(r).^2/2/pi;
    
end

uR = [real(uR); imag(uR)]


end