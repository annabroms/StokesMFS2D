function selfTestCapInterp()
%SELFTESTCAPINTERP Lightweight self-test of the interpolation machinery used
% by sep25_capacitance_cmap_gap_interp.m. Uses only synthetic, analytically
% known matrix functions, so it runs in well under a second and does not
% touch the Laplace MFS code at all: it checks the generic numerical
% machinery (Chebyshev nesting, barycentric interpolation accuracy and
% exactness at nodes, common-basis reconstruction, numerical rank), not the
% physics of the pair map itself.
%
% Anna Broms, Sep 2026

fprintf('=== selfTestCapInterp ===\n');
n_fail = 0;

%% 1. Nesting of Chebyshev-Lobatto grids for q=9,17,33,65
smin = -3; smax = 1;
q_levels = [9 17 33 65];
s_prev = [];
for q = q_levels
    s = chebLobattoNodes(smin,smax,q);
    assert(abs(s(1)-smax)<1e-12 && abs(s(end)-smin)<1e-12, ...
        'chebLobattoNodes: endpoints wrong.');
    if ~isempty(s_prev)
        d = min(abs(s(:)-s_prev(:)'),[],1); % nearest new-grid node to each old node
        ok = all(d < 1e-9*(smax-smin));
        n_fail = n_fail + report(sprintf('nesting q=%d contains previous grid',q),~ok,q);
    end
    s_prev = s;
end

%% 2. Barycentric interpolation reproduces a smooth analytic matrix function
% F(s) = exp(s)*M1 + sin(2*s)*M2, a smooth (entire) function of s, so a
% degree-64 Chebyshev interpolant should be accurate to near machine
% precision, and interpolation error should shrink as q grows.
rng(0);
M1 = randn(5,4); M2 = randn(5,4);
Ffun = @(s) exp(s)*M1 + sin(2*s)*M2;

s_test = linspace(smin,smax,37)'; % includes non-node points
err_q = nan(size(q_levels));
for i = 1:numel(q_levels)
    q = q_levels(i);
    s_nodes = chebLobattoNodes(smin,smax,q);
    w = chebBarycentricWeights(q);
    F_nodes = zeros(5,4,q);
    for j = 1:q
        F_nodes(:,:,j) = Ffun(s_nodes(j));
    end
    Fq = evalMatrixChebBary(s_test,s_nodes,w,F_nodes);
    e = 0;
    for k = 1:numel(s_test)
        e = max(e,norm(Fq(:,:,k)-Ffun(s_test(k)),'fro'));
    end
    err_q(i) = e;
end
n_fail = n_fail + report('interpolation error decreases with q',~all(diff(err_q)<=1e-3*err_q(1:end-1)+1e-13),err_q);
n_fail = n_fail + report('finest-level interpolation is accurate (< 1e-8)',err_q(end)>1e-8,err_q(end));

% Exactness at a training node.
q = q_levels(end);
s_nodes = chebLobattoNodes(smin,smax,q);
w = chebBarycentricWeights(q);
F_nodes = zeros(5,4,q);
for j = 1:q
    F_nodes(:,:,j) = Ffun(s_nodes(j));
end
Fexact = evalMatrixChebBary(s_nodes(3),s_nodes,w,F_nodes);
n_fail = n_fail + report('exact reproduction at a training node', ...
    norm(Fexact-F_nodes(:,:,3),'fro')>1e-10*max(1,norm(F_nodes(:,:,3),'fro')),NaN);

%% 3. Common-basis reconstruction is exact when r_max = full rank
n = 6; q = 5;
dA = zeros(n,n,q);
for j = 1:q
    dA(:,:,j) = randn(n)*sin(j); % rank-n generic snapshots
end
basis = commonBasisFromSnapshots(dA,n);
max_err = 0;
for j = 1:q
    P = basis.U*(basis.U'*dA(:,:,j)*basis.V)*basis.V';
    max_err = max(max_err,norm(P-dA(:,:,j),'fro')/max(1,norm(dA(:,:,j),'fro')));
end
n_fail = n_fail + report('full-rank common basis reconstructs snapshots exactly', ...
    max_err>1e-8,max_err);

%% 4. numericalRank basic sanity
sigma = [10 5 1 1e-6 1e-12]';
r = numericalRank(sigma,[1e-3 1e-9]);
n_fail = n_fail + report('numericalRank([...],[1e-3 1e-9]) == [3 4]', ...
    ~isequal(r,[3 4]),r);

if n_fail == 0
    fprintf('selfTestCapInterp: ALL CHECKS PASSED.\n\n');
else
    error('selfTestCapInterp:Failed','%d check(s) failed, see above.',n_fail);
end

end

function n_fail = report(name,failed,val)
if failed
    fprintf('  [FAIL] %s (value: %s)\n',name,mat2str(val));
    n_fail = 1;
else
    fprintf('  [ OK ] %s\n',name);
    n_fail = 0;
end
end
