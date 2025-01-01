function  [Y,Bi1]  = getSVDFactors(N,tol,visualise)

%econ = 0; 

% if econ
%     [UU,S,V] = svd(N,"econ");
%     S = diag(S);
%     iS = 1./S; % rank
%     Y = V*diag(iS); 
%     Bi1 = UU; 
%else

[UU,S,V] = svd(N);
S = diag(S);

%use relative tolerance 
ra = sum(S>max(S)*tol); 

%with tol = eps, we throw away only the last singval, corresponding to the
%non-trivial null-space (normal direction for single layer...). Does that make sense?

%yes! testing
%ra = size(N,2); %leads to a horribly irregular 1-body basis. 


if visualise
   SS = diag(S);
   % Sv = vecnorm(SS*V',2,1);

    figure(57)
   % clf;
    semilogy(S,'o-');
    hold on
    semilogy(ra*ones(1,2),logspace(-15,5,2),'r--')
    c = max(S)/min(S);  %Condition number
    str = sprintf('Self condition number %1.3e, max sing %1.3e',c,max(S));
    title(str,'interpreter','latex');
    grid on
    xlabel('$j$','interpreter','latex')
    ylabel('$\sigma_j$','interpreter','latex')

    % figure(56)
    % SV = SS*V';
    % semilogy(abs(SV(:,end-50:end)'))
    % 
    % 
    % figure(56)
    % semilogy(abs(UU(end-3:end,:)'))

    % figure(57)
    % semilogy(abs(diff(S)));
    % hold on
    % title('Decay rate of sing vals','interpreter','latex')
    % 
    % 
    % figure(58)
    % semilogy(abs(diff(S)./S(2:end)));
    % hold on
    % title('Relative decay rate of sing vals','interpreter','latex')
end

S = S(1:ra);  
iS = 1./S; % rank
Y = V(:,1:ra)*diag(iS); 
Bi1 = UU(:,1:ra); 



%invert matrix containing singular values, up to a tolerance





end