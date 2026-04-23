% Compare looped and batched one-body transforms in Laplace peanut matvecs.
%
% The batched path preserves the two-stage ordering U{1} first, Y{1}
% second.  For elastance/projected solves, the per-particle charge projection 
% is applied columnwise after the batched one-body map.

script_name = mfilename;
script_date = 'Apr 23, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (%s) ===\n',script_name,script_date);

if ~exist('P_sweep','var') || isempty(P_sweep)
    P_sweep = [10 40 100 250 1000];
end
if ~exist('R','var') || isempty(R)
    R = 2;
end
if ~exist('N_c','var') || isempty(N_c)
    N_c = 40;
end
if ~exist('N_f','var') || isempty(N_f)
    N_f = 80;
end
if ~exist('n_repeats','var') || isempty(n_repeats)
    n_repeats = 3;
end
if ~exist('rng_seed','var') || isempty(rng_seed)
    rng_seed = 230423;
end

opt = getLaplace2Dparams(max(P_sweep),R,N_c,N_f);
a_c = opt.a_c;

tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1).';
rbase_out_c = R*(cos(tout_c) + 1i*sin(tout_c));
N_large = numel(rbase_out_c);

tin_c = linspace(0,2*pi,N_c+1).';
tin_c = tin_c(1:end-1);
rbase_in_c = opt.Rp_c*(cos(tin_c) + 1i*sin(tin_c));

project_modes = [false true];
rng(rng_seed,'twister');

for imode = 1:numel(project_modes)
    project_charge = project_modes(imode);
    [U,Y] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c, ...
        [0 N_large],project_charge);
    Nii = lapSLPmat(rbase_in_c,rbase_out_c);
    self_results = zeros(numel(P_sweep),5);

    fprintf('\nproject_charge=%d  P sweep [%s], N_c=%d, N_check=%d, repeats=%d\n', ...
        project_charge,num2str(P_sweep),N_c,N_large,n_repeats);
    fprintf('One-body transform: lam = Y{1}*(U{1}*tau)\n');
    fprintf('%8s %12s %12s %10s %12s\n', ...
        'P','loop_s','batched_s','speedup','rel_inf');

    for ip = 1:numel(P_sweep)
        P = P_sweep(ip);
        tau = randn(P*N_large,1);

        loop_result = apply_one_body_loop( ...
            tau,U,Y,P,N_large,N_c,project_charge);
        batch_result = apply_one_body_batched( ...
            tau,U,Y,P,N_large,N_c,project_charge);
        rel_inf = norm(batch_result-loop_result,inf) / ...
            max(1,norm(loop_result,inf));

        t_loop = zeros(n_repeats,1);
        t_batch = zeros(n_repeats,1);
        for rep = 1:n_repeats
            t_loop(rep) = timeit(@() apply_one_body_loop( ...
                tau,U,Y,P,N_large,N_c,project_charge));
            t_batch(rep) = timeit(@() apply_one_body_batched( ...
                tau,U,Y,P,N_large,N_c,project_charge));
        end

        loop_s = median(t_loop);
        batch_s = median(t_batch);
        speedup = loop_s / max(batch_s,eps);
        fprintf('%8d %12.4e %12.4e %10.2f %12.3e\n', ...
            P,loop_s,batch_s,speedup,rel_inf);

        res0 = randn(P*N_large,1);
        self_loop_result = apply_self_block_loop(res0,loop_result,Nii,P,N_large,N_c);
        self_batch_result = apply_self_block_batched(res0,loop_result,Nii,P,N_large,N_c);
        self_rel_inf = norm(self_batch_result-self_loop_result,inf) / ...
            max(1,norm(self_loop_result,inf));

        t_self_loop = zeros(n_repeats,1);
        t_self_batch = zeros(n_repeats,1);
        for rep = 1:n_repeats
            t_self_loop(rep) = timeit(@() apply_self_block_loop( ...
                res0,loop_result,Nii,P,N_large,N_c));
            t_self_batch(rep) = timeit(@() apply_self_block_batched( ...
                res0,loop_result,Nii,P,N_large,N_c));
        end

        self_loop_s = median(t_self_loop);
        self_batch_s = median(t_self_batch);
        self_speedup = self_loop_s / max(self_batch_s,eps);
        self_results(ip,:) = [P,self_loop_s,self_batch_s, ...
            self_speedup,self_rel_inf];
    end

    fprintf('Self-block correction: res = res - Nii*lam_self per body\n');
    fprintf('%8s %12s %12s %10s %12s\n', ...
        'P','loop_s','batched_s','speedup','rel_inf');
    for ip = 1:size(self_results,1)
        fprintf('%8d %12.4e %12.4e %10.2f %12.3e\n', ...
            self_results(ip,1),self_results(ip,2),self_results(ip,3), ...
            self_results(ip,4),self_results(ip,5));
    end
end

function lambda = apply_one_body_loop(tau,U,Y,P,N_large,N_c,project_charge)
lam_c = zeros(P*N_c,1);

for i = 1:P
    block = (i-1)*N_large+1:i*N_large;
    tau_i = tau(block);
    lam_i_nonp = Y{1}*(U{1}*tau_i);
    lam_i = project_charge_mode(lam_i_nonp,project_charge);
    idx = (i-1)*N_c+1:i*N_c;
    lam_c(idx) = lam_i;
end

lambda = lam_c;
end

function lambda = apply_one_body_batched(tau,U,Y,P,N_large,N_c,project_charge)
tau_blocks = reshape(tau(1:P*N_large),N_large,P);
lam_blocks_nonp = Y{1}*(U{1}*tau_blocks);
if project_charge
    lam_blocks = lam_blocks_nonp - mean(lam_blocks_nonp,1);
else
    lam_blocks = lam_blocks_nonp;
end
lambda = reshape(lam_blocks(1:N_c,:),[],1);
end

function res = apply_self_block_loop(res,lam_self,Nii,P,N_large,N_c)
for i = 1:P
    idx = (i-1)*N_c+1:i*N_c;
    block = (i-1)*N_large+1:i*N_large;
    uii = Nii*lam_self(idx);
    res(block) = res(block) - uii;
end
end

function res = apply_self_block_batched(res,lam_self,Nii,P,N_large,N_c)
lam_self_blocks = reshape(lam_self,N_c,P);
uii_blocks = Nii*lam_self_blocks;
res = res - reshape(uii_blocks(1:N_large,:),[],1);
end

function lam_out = project_charge_mode(lam_in,project_charge)
if ~project_charge || isempty(lam_in)
    lam_out = lam_in;
    return
end
n = numel(lam_in);
lam_out = lam_in - (sum(lam_in)/n);
end
