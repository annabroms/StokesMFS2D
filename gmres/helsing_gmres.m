function [x,it,resvec,trueres] = helsing_gmres(f,b,n,m,tol,opt,grid)
% *** GMRES with low-threshold stagnation control taken from RCIP code of Helsing***
% Optional verbosity input is passed in opt.gmres_verbose:
%   0 = silent, 1 = final summary, 2 = per-iteration + final summary
%
% Backward compatibility:
%   If the 6th input is numeric, it is interpreted as gmres verbosity.
if nargin < 6
    opt = struct();
end
if nargin < 7
    grid = [];
end

if isnumeric(opt) && isscalar(opt)
    gmres_verbose = opt;
elseif isstruct(opt) && isfield(opt,'gmres_verbose') && ...
        ~isempty(opt.gmres_verbose) && isnumeric(opt.gmres_verbose) && ...
        isscalar(opt.gmres_verbose)
    gmres_verbose = opt.gmres_verbose;
else
    gmres_verbose = 0;
end

% Debug levels:
%   0: silent
%   1: final summary only
%   2: per-iteration residuals + final summary
verbose_final = gmres_verbose >= 1;
verbose_iter = gmres_verbose >= 2;

V=zeros(n,m+1);
H=zeros(m);
cs=zeros(m,1);
sn=zeros(m,1);
bnrm2=norm(b);
if bnrm2 == 0
    x = zeros(n,1);
    it = 0;
    resvec = 0;
    trueres = 0;
    if verbose_final
        fprintf('GMRES summary: gmres_iter=%d, est_relres=%.3e, true_relres=%.3e\n', ...
            it, 0.0, 0.0);
    end
    return
end

V(:,1)=b/bnrm2;
s=bnrm2*eye(m+1,1);
for it=1:m  
it1=it+1;                                   
%w=A*V(:,it);
w=f(V(:,it))-V(:,it);
for k=1:it
  H(k,it)=V(:,k)'*w;
  w=w-H(k,it)*V(:,k);
end
H(it,it)=H(it,it)+1;
wnrm2=norm(w);
V(:,it1)=w/wnrm2;
for k=1:it-1                                
  temp     = cs(k)*H(k,it)+sn(k)*H(k+1,it);
  H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
  H(k,it)  = temp;
end
[cs(it),sn(it)]=rotmat(H(it,it),wnrm2);     
H(it,it)= cs(it)*H(it,it)+sn(it)*wnrm2;
s(it1) =-sn(it)*s(it);                      
s(it)  = cs(it)*s(it);                         
est_relres = abs(s(it1))/bnrm2;
resvec(it) = est_relres;
    if verbose_iter
        fprintf('GMRES iter %4d: est_relres = %.3e\n',it,est_relres);
    end





if (est_relres<=tol)||(it==m)                     
  y=triu(H(1:it,1:it))\s(1:it);             
  x=fliplr(V(:,1:it))*flipud(y);
  %trueres=norm(A*x-b)/bnrm2;
  trueres=norm(f(x)-b)/bnrm2;
  if verbose_final
      fprintf('GMRES summary: gmres_iter=%d, est_relres=%.3e, true_relres=%.3e\n', ...
          it, est_relres, trueres);
  end
  resvec = resvec(1:it);
  break
end
end
end
