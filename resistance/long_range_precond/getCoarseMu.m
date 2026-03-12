function mu_coarse = getCoarseMu(f,Sinv,Z,Y,db,P,Nc,a,Vk)

if nargin<9

    proj_rhs = zeros(db*P,1);
    for k = 1:P
       % proj_rhs(db*(k-1)+1:db*(k-1)+db/2) = Mx'*f((k-1)*Nc*a+1:k*Nc*a); %x part
       % proj_rhs(db*(k-1)+db/2+1:db*k) = Mx'*f(P*Nc*a+(k-1)*Nc*a+1:P*Nc*a+k*Nc*a); %y part
    
        proj_rhs(db*(k-1)+1:db*k) = Y'*[f((k-1)*Nc*a+1:k*Nc*a); f(P*Nc*a+(k-1)*Nc*a+1:P*Nc*a+k*Nc*a)]; %x part, y part
    end
    coarse_coeff = Sinv*proj_rhs; % This matvec can be large if many basis functions
    mu_coarse = zeros(2*Nc*P*a,1); 
    for k = 1:P
       % tau_coarse((k-1)*Nc+1:k*Nc) = Nx*coarse_coeff((k-1)*db+1:(k-1)*db+db/2);
       % tau_coarse((k-1)*Nc+1+P*Nc:k*Nc+P*Nc) = Ny*coarse_coeff((k-1)*db+db/2+1:k*db);
    
       res = Z*coarse_coeff((k-1)*db+1:k*db);
       mu_coarse((k-1)*Nc*a+1:k*Nc*a) = res(1:end/2);
       mu_coarse((k-1)*Nc*a+1+P*Nc*a:k*Nc*a+P*Nc*a) = res(end/2+1:end);
    
    end
else
    mu_coarse = Vk*Sinv*Vk'*f;
end

end