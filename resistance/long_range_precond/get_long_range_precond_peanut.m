function [Sinv,Zi,Yi,db] = get_long_range_precond_peanut(q,rin,rout,rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,opt)
%%GET_LONG_RANGE_PRECOND_MU  Construct coarse-space projection matrices for long-range preconditioning.
% The same thing as before, but for mu instead of lambda.
%
%   [Sinv, Zi, Yi,db] = GET_LONG_RANGE_PRECOND_PEANUT(q, rin, rout, opt)
%
%   Constructs the matrices used in the long-range preconditioner:
%   the coarse-to-fine mappings AN and AM, and the inverse coarse interaction matrix Rinv.
%   These are used to define the projectors:
%
%       P = I - G * Z * Sinv * Y'     (left-side projection, see applyPmat)
%       Q = I - Z * Sinv * Y' * G     (right-side projection, see applyQmat)
%
%   The choice of coarse basis (e.g. translation-only vs. rigid-body modes) 
%   is set via the option opt.lr.
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
%     - The function assumes block structure per body and Yi and Zi constructs
%       block-diagonal matrices Z and Y accordingly.
%
%   See also: applyPmat, applyQmat, get_long_range_precond_mob
%
% Anna Broms, Oct 19, 2025

lr = opt.lr;
a = opt.a_c;
Nc = opt.N_c; 
P = length(q);  


%% Is it a good choice to take singular vectors of G for a single body to be the coarse basis functions in Z, Y? Yes!

%svd_type

G = stokSLPmat(rin(1:Nc),rout(1:Nc*a),1); %rectangular matrix 
[Uii,Yii] = getSelfPseudo(1,rin(1:Nc),rout(1:Nc*a),[],[],[1,opt.N_c*opt.a_c],opt.s);
%[V,D] = eig(G);
[U,S,V] = svd(G); 
smax = lr-2;
%s = diag(S);
%Zi = V(:,1:smax);%*diag(1./s(1:smax));
%Zi = Yii{1}*(Uii{1}*U(:,1:smax));
Zi = U(:,1:smax);

%V = [ones(Nc,1) zeros(Nc,1); zeros(Nc,1) ones(Nc,1)];
%smax = 2;
Zi = G*V(:,1:smax);
%Zi = [ones(Nc*a,1) zeros(Nc*a,1); zeros(Nc*a,1) ones(Nc*a,1)];
Yi = U(:,1:smax);
%Yi = Zi; %testing
db = smax;



%% call fmm here instead with 2PK sources per body, all to all
%put this in an eval multiD function instead!

targ = [real(rout)';imag(rout)'];



%% Build matrix Vmat using its block structure instead. 

% Exteremely slow implementation: 
S = zeros(db*P);
Vmat = zeros(Nc*P*2*a,P*db);
mu = zeros(Nc*P*2*a,1);
for k = 1:P
    k
    mu(:) = 0;
    
    
    for i = 1:db
        mu((k-1)*Nc*a+1:k*Nc*a) = Zi(1:end/2,i);
        mu((k-1)*Nc*a+1+Nc*P*a:k*Nc*a+Nc*P*a) = Zi(end/2+1:end,i);

        res = matvec_res_peanut(mu,rbase_in_c,rbase_in_f,rin,rbase_out_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,0);
        Vmat(:,db*(k-1)+i) = res;
    end

    
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


%Make a coarse version of the same thing
% Gc = stokSLPmat(rin_c,rout_c,1);
% if lr == 2
%     AMc = getGlobalKmat2D(rout_c,q);
%     ANc= getGlobalKmat2D(rin_c,q);
% else
%     ANc = [blkdiag(blocksNx_c{:}); blkdiag(blocksNy_c{:})]./Nc_c;
%     AMc = [blkdiag(blocksMx_c{:}); blkdiag(blocksMy_c{:})]./Nc_c/a;
% end
% 
% Rc = AMc'*Gc*ANc;
% Rinv = (Rc\eye(db*P)); %This might be off by some scaling. 


%B = (R\eye(2*P))*AM';

%B = Rinv*AM';
%B = (Rc\eye(3*P))*AM';%us coarse R instead
%Cmap = AN*B;
%Pmat = eye(round(2*P*a*Nc))-G*Cmap;
%Qmat = eye(round(2*P*Nc))-Cmap*G;

%Convergence for the system involving R?
% e = eig(R);
% figure()
% plot(real(e),imag(e),'+')
% 
% ob_block = Y'*G*Z; %one_body_block;
% ob_inv = ob_block\eye(db);
% B = kron(eye(P), ob_inv);
% R2 = R*B;
% e = eig(R2);
% hold on
% plot(real(e),imag(e),'ro')

%Zi = G*Zi;

%Zi = U(:,1:smax); 
 
end

