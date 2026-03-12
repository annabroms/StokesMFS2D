function [Sinv,Zi,Yi,db] = get_long_range_precond_mob(q,rin,rout,L,Lr,opt)
%%GET_LONG_RANGE_PRECOND_MOB  Construct coarse-space projection matrices for long-range preconditioning.
%
%   [Sinv, Z, Y,db] = GET_LONG_RANGE_PRECOND_MOB(q, rin, rout, L,Lr,opt)
%
%   Constructs the matrices used in the long-range preconditioner:
%   the diagonal blocks in the matrices Y and Z, and the inverse coarse interaction matrix Sinv.
%   These are used to define the projectors:
%
%       P = I - G_L * Z * Sinv * Y'     (left-side projection, see applyPmat_mob)
%       Q = I - Zi * Sinv * Yi' * G_L     (right-side projection, see applyQmat_mob)
%
%   The choice of coarse basis is set via the option opt.lr.
%
%   INPUTS:
%     q        - P array of particle centers.
%     rin      - PN complex valued array of proxy/source points (used here to determine for coarse basis flow).
%     rout     - PM complex valued array of surface/target points (used to evaluate coarse flow).
%     opt      - Struct with fields:
%                  • lr: 1 = translation-only modes (2 per body),
%                        2 = full rigid-body modes (3 per body),
%                        >3 use columns of the single particle SVD
%                  • N: number of proxy points per body,
%                  • M: number of collocation points per body,
%                  • other kernel options (such as fmm flat passed to getFlow).
%
%   OUTPUTS:
%     Sinv     - (Pk×Pk) inverse of coarse interaction matrix S = blkdiag(Y')*U,
%                where U is the coarse-flow matrix at the PM targets from
%                the Pk coarse basis functions
%     Z       - (2PM×Pk) matrix mapping coarse coefficients to proxy source strengths.
%     Y       - (2PN×Pk) matrix whose transpose, Y', maps to coarse velocity space.
%
%   NOTES:
%     - The function assumes block structure per body 
%
%   See also: applyPmat_mob, applyQmat_mob
%
% Anna Broms, Oct 22, 2025

lr = opt.lr;
a = opt.a_c;
Nc = opt.N_c; 
P = length(q);  

%Use different matrices Z and Y using ones matrices instead!
if lr == 1
    for k = 1:P
        blocksNx{k} = [ones(Nc,1) zeros(Nc,1)];
        blocksNy{k} = [zeros(Nc,1) ones(Nc,1)];
        blocksMx{k} = [ones(Nc*a,1) zeros(Nc*a,1)];
        blocksMy{k} = [zeros(Nc*a,1) ones(Nc*a,1)];
        
        % %Coarse version of the same thing
        % blocksNx_c{k} = [ones(Nc_c,1) zeros(Nc_c,1)];
        % blocksNy_c{k} = [zeros(Nc_c,1) ones(Nc_c,1)];
        % blocksMx_c{k} = [ones(Nc_c*a,1) zeros(Nc_c*a,1)];
        % blocksMy_c{k} = [zeros(Nc_c*a,1) ones(Nc_c*a,1)];
       
    end
    %AN = [blkdiag(blocksNx{:}); blkdiag(blocksNy{:})]./Nc; %block diagonal matrix with all the Y_i
    %AM = [blkdiag(blocksMx{:}); blkdiag(blocksMy{:})]./Nc/a;%block diagonal matrix with all the Z_i
    db = 2; %dimension per body of the coarse source
    Kout = [ones(Nc*a,1) zeros(Nc*a,1); zeros(Nc*a,1) ones(Nc*a,1)]./Nc/a;
    svd_type = 0; 

elseif lr == 2
    %AM = getGlobalKmat2D(rout,q); %block diagonal matrix with all the Y_i
    %AN = getGlobalKmat2D(rin,q);   %block diagonal matrix with all the Z_i
    db = 3; %dimension per body of the coarse source
    Kout = getKmat2D(rout(1:Nc*a),q(1));
    svd_type = 0; 
elseif lr>2 
    lmax = lr-3;
    %lmax = 0; 
    
    % Start with first two columns: (1,0) for Nx, (0,1) for Ny
    % Nx = [ones(Nc,1), zeros(Nc,1)];
    % Ny = [zeros(Nc,1), ones(Nc,1)];
    % Mx = [ones(Nc*a,1), zeros(Nc*a,1)];
    % My = [zeros(Nc*a,1), ones(Nc*a,1)];
    svd_type = 1; 
    if ~svd_type
        %use lmax fourier modes. lmax =0 should give back the results from lr =1.
        tN = linspace(0,2*pi,Nc+1);
        tN = tN(1:end-1)'; 
        tM = linspace(0,2*pi,Nc*a+1);
        tM = tM(1:end-1)';
        Nx = [ones(Nc,1)];
      %  Ny = [ones(Nc,1)];
        Mx = [ones(Nc*a,1)];
       % My = [ones(Nc*a,1)];
        
        if lmax >0
            Nx = [Nx, sin(tN)]; % will use the same basis in the x and y component
            Mx = [Mx, sin(tM)];
    
            % Add harmonic blocks up to order lmax
            for l = 2:lmax
                beta = l^2;
                %beta = 1; 
                Nx = [Nx, beta*sin((l)*tN), beta*cos((l)*tN)];
                Mx = [Mx, sin((l)*tM), cos((l)*tM)];         
            end
        end
    
        
    
    
        for k = 1:P
            % Store
            blocksNx{k} = [Nx zeros(size(Nx))];
            blocksNy{k} = [zeros(size(Nx)) Nx];
          %  blocksMx{k} = Mx;
            %blocksMy{k} = My;
            blocksMx{k} = [Mx zeros(size(Mx))];
            blocksMy{k} = [zeros(size(Mx)) Mx];
            
        end
    
        
    
        if lmax
            db = (lmax-1)*4+4; 
        else
            db = 2;
        end
       % AN = [blkdiag(blocksNx{:}); blkdiag(blocksNy{:})];
        %AM = [blkdiag(blocksMx{:}); blkdiag(blocksMy{:})];
        Ny = Nx;
    
    end
    
end

%% 
if svd_type
    G = singleLayer(rin(1:Nc),rout(1:Nc*a),1); %rectangular matrix 
    %[V,D] = eig(G);
    %[U,S,V] = svd(G*(eye(2*Nc)-L)+Lr); 
    [U,S,V] = svd(G*(eye(2*Nc)-L)); 
    %[U,S,V] = svd(G); 
    smax = lr-2;
    s = diag(S);
    Zi = V(:,1:smax);%*diag(1./s(1:smax));
    Yi = U(:,1:smax);
    Nx = []; Ny = []; Mx = []; 
    
    % figure()
    % d = diag(D); 
    % plot(real(d),imag(d),'*');
    % %R1 = AM'*G*AN;
    % max_eig = 1;
    % [B,I] = maxk(abs(d),max_eig);
    % Vk = V(:,I);
    % figure()
    % subplot(1,2,1)
    % plot(real(Vk));
    % subplot(1,2,2);
    % plot(imag(Vk));
    % Nxx = []; Nyy = []; Nxy = []; Nyx = [];
    % for k = 1:max_eig
    %     Nxx = [Nxx real(Vk(1:end/2,k))];
    %     %if norm(imag(Vk(1:end/2,k)))>1e-4
    %         Nxy = [Nxy imag(Vk(1:end/2,k))];
    %     %end
    %     % else
    %     %     disp('here x');
    %     % end
    %     Nyy = [Nyy real(Vk(end/2+1:end,k))];
    %     %if norm(imag(Vk(end/2+1:end,k)))>1e-4
    %         Nyx = [Nyx imag(Vk(end/2+1:end,k))];
    %     %end
    %     %     disp('here y');
    %     % end
    % end
    % 
    % NN = [Nxx Nxy; Nyx Nyy];
    
    db = smax;
% 
% Mx = Nx; %only works if the same number of source and target points
% Create som proxy for Mx
% Mx = [ones(Nc*a,1)]; lmax = db;  
% if lmax >0
%     Nx = [Nx, sin(tN)]; % will use the same basis in the x and y component
%     Mx = [Mx, sin(tM)];
% 
%     % Add harmonic blocks up to order lmax
%     for l = 2:lmax
%         beta = l^2;
%         %beta = 1; 
%         Nx = [Nx, beta*sin((l)*tN), beta*cos((l)*tN)];
%         Mx = [Mx, sin((l)*tM), cos((l)*tM)];         
%     end
% end

%Yp = [Mx zeros(size(Mx)); zeros(size(Mx)) Mx];
%Y = MM; 
else
    Yi = [Mx zeros(size(Mx)); zeros(size(Mx)) Mx];
    Zi = [Nx zeros(size(Mx)); zeros(size(Mx)) Ny];
end

%% call fmm here instead with 2PK sources per body, all to all
%put this in an eval multiD function instead!
eps = 1e-10; %must be chosen sufficiently small!
targ = [real(rout)';imag(rout)'];
% srcinfo.sources = [real(rin)'; imag(rin)'];    
% src_struct = zeros(db*P,2,P*Nc);
% for k = 1:db*P
%     src_struct(k,1,:) = AN(1:end/2,k)';
%     src_struct(k,2,:) = AN(end/2+1:end,k)';
% end
% 
% srcinfo.stoklet = src_struct; %should be nd,2,n
% srcinfo.nd = db*P; 
% Ufmm = stfmm2d(eps, srcinfo, 0, targ, 1); % might want to change eps here
% Vmat2 = [reshape(Ufmm.pottarg(:,1,:),srcinfo.nd,Nc*P*a) reshape(Ufmm.pottarg(:,2,:),srcinfo.nd,Nc*P*a)]'./2/pi;
%R = AM'*Vmat; %This can certainly be computed faster! 

%% Build matrix Vmat using its block structure instead. 
% Enough to use sources on one body at the time due to the sparsity structure
S = zeros(db*P);
Vmat = zeros(Nc*P*2*a,P*db);
LrZ = Lr*Zi;
LZ = Zi-L*Zi;
LZ = Zi;
% Let the coarse operator be the pure SLP or G(I-L)+Lr?

for k = 1:P
    %Do a block call for the 2k basis functions on this body
    srcinfo.sources = [real(rin((k-1)*Nc+1:k*Nc))'; imag(rin((k-1)*Nc+1:k*Nc))'];   
    srcinfo.nd = db;
    src_struct = zeros(db,2,Nc);

    for i = 1:db
        if svd_type
            
            src_struct(i,1,:) = LZ(1:end/2,i);
            src_struct(i,2,:) = LZ(end/2+1:end,i);
        else
            % Assuming basis functions with no coupling between x and y,
            % such as the constand basis, or using strategy 3, with
            % trigonometric modes. 
            if i < db/2 +1
                src_struct(i,1,:) = Nx(:,i)';
                src_struct(i,2,:) = zeros(Nc,1);
            else
                src_struct(i,1,:) = zeros(Nc,1);
                src_struct(i,2,:) = Ny(:,i-db/2)';
            end
        end
   
    end

    srcinfo.stoklet = src_struct; %should be nd,2,n
    srcinfo.nd = db; 
    Ufmm = stfmm2d(eps, srcinfo, 0, targ, 1); 
    if smax>1
        x_comp = reshape(Ufmm.pottarg(:,1,:),srcinfo.nd,Nc*P*a)'./2/pi;
        y_comp = reshape(Ufmm.pottarg(:,2,:),srcinfo.nd,Nc*P*a)'./2/pi;
    else
        x_comp = Ufmm.pottarg(1,:);
        y_comp = Ufmm.pottarg(2,:);
    end    

    %replace with direct evaluation instead
   % [x_comp,y_comp] = stokesletDirect(real(rin((k-1)*Nc+1:k*Nc))',imag(rin((k-1)*Nc+1:k*Nc))',...
    %    real(rout),imag(rout),Z(1:end/2,:),Z(end/2+1:end,:),Nc);
   % x_comp(Nc*a*(k-1)+1:Nc*a*k,:) = x_comp(Nc*a*(k-1)+1:Nc*a*k,:)+LrZ(1:end/2,:);
    %y_comp(Nc*a*(k-1)+1:Nc*a*k,:) = y_comp(Nc*a*(k-1)+1:Nc*a*k,:)+LrZ(end/2+1:end,:);
    Vmat(1:end/2,db*(k-1)+1:db*k) = x_comp;
    Vmat(end/2+1:end,db*(k-1)+1:db*k) = y_comp;

   %Slow, but requires less memory. 
    % for i = 1:P
    %     R((i-1)*db+1:i*db,(k-1)*db+1:k*db) = R((i-1)*db+1:i*db,(k-1)*db+1:k*db) + ...
    %         Y'*[x_comp((i-1)*Nc*a+1:i*Nc*a,:); y_comp((i-1)*Nc*a+1:i*Nc*a,:)];
    % end
    
end

% Check "loss of information" in Vmat
% [U,S,V] = svd(Vmat);
% ss = diag(S);
% figure(33)
% semilogy(ss)
% hold on

for k = 1:P
    S((k-1)*db+1:k*db,:) = Yi'*[Vmat((k-1)*Nc*a+1:k*Nc*a,:); Vmat(P*Nc*a+(k-1)*Nc*a+1:P*Nc*a+k*Nc*a,:)]; 
end

disp('... Matrix S is built')

Sinv = (S\eye(db*P));


    
end

