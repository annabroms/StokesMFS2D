% Compare looped and batched one-body work in resistance matvecs.
%
% The batched path preserves the two-stage ordering U{1} first, Y{1}
% second. It deliberately does not form Y{1}*U{1}.

script_name = mfilename;
script_date = 'Apr 23, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (%s) ===\n',script_name,script_date);

P_sweep = [10 40 100 250 1000];
rng_seed = 230423;
N_c = 60;
n_repeats = 10;

opt = get2Dparams(max(P_sweep),N_c,N_f);
a_c = opt.a_c;

tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1).';
rbase_out_c = cos(tout_c) + 1i*sin(tout_c);
N_large = numel(rbase_out_c);

tin_c = linspace(0,2*pi,N_c+1).';
tin_c = tin_c(1:end-1);
rbase_in_c = opt.Rp_c*cos(tin_c) + 1i*opt.Rp_c*sin(tin_c);

svd_opts = struct( ...
    'column_weight',logical(getOptField(opt,'column_weight',false)), ...
    'left_weight',logical(getOptField(opt,'left_weight',false)));
[U,Y] = getSelfPseudo(1,rbase_in_c,rbase_out_c,[],[], ...
    [0,N_large],0,svd_opts);
Nii = stokSLPmat(rbase_in_c,rbase_out_c,1);

fprintf(['One-body sizes: P sweep [%s], N_c=%d, N_check=%d, ', ...
    'U=%dx%d, Y=%dx%d, repeats=%d\n'], ...
    num2str(P_sweep),N_c,N_large, ...
    size(U{1},1),size(U{1},2),size(Y{1},1),size(Y{1},2), ...
    n_repeats);
fprintf('\nOne-body transform: lambda = Y{1}*(U{1}*tau)\n');
fprintf('%8s %12s %12s %10s %12s\n', ...
    'P','loop_s','batched_s','speedup','rel_inf');

rng(rng_seed,'twister');
transform_results = repmat(struct('P',0,'loop_s',0,'batch_s',0, ...
    'speedup',0,'rel_inf',0),numel(P_sweep),1);
self_results = transform_results;

for ip = 1:numel(P_sweep)
    P = P_sweep(ip);
    tau = randn(2*P*N_large,1);

    loop_result = apply_one_body_loop(tau,U,Y,P,N_large,N_c);
    batch_result = apply_one_body_batched(tau,U,Y,P,N_large,N_c);
    rel_inf = norm(batch_result-loop_result,inf) / ...
        max(1,norm(loop_result,inf));

    t_loop = zeros(n_repeats,1);
    t_batch = zeros(n_repeats,1);
    for rep = 1:n_repeats
        t_loop(rep) = timeit(@() apply_one_body_loop( ...
            tau,U,Y,P,N_large,N_c));
        t_batch(rep) = timeit(@() apply_one_body_batched( ...
            tau,U,Y,P,N_large,N_c));
    end

    loop_s = median(t_loop);
    batch_s = median(t_batch);
    speedup = loop_s / max(batch_s,eps);
    fprintf('%8d %12.4e %12.4e %10.2f %12.3e\n', ...
        P,loop_s,batch_s,speedup,rel_inf);

    transform_results(ip) = struct('P',P,'loop_s',loop_s, ...
        'batch_s',batch_s,'speedup',speedup,'rel_inf',rel_inf);

    PM = P*N_large;
    res0 = randn(2*PM,1);
    self_loop_result = apply_self_block_loop( ...
        res0,loop_result,Nii,P,N_large,N_c);
    self_batch_result = apply_self_block_batched( ...
        res0,loop_result,Nii,P,N_large,N_c);
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
    self_results(ip) = struct('P',P,'loop_s',self_loop_s, ...
        'batch_s',self_batch_s,'speedup',self_speedup, ...
        'rel_inf',self_rel_inf);
end

fprintf('\nSelf-block subtraction: res = res - Nii*lambda_self per body\n');
fprintf('%8s %12s %12s %10s %12s\n', ...
    'P','loop_s','batched_s','speedup','rel_inf');
for ip = 1:numel(self_results)
    fprintf('%8d %12.4e %12.4e %10.2f %12.3e\n', ...
        self_results(ip).P,self_results(ip).loop_s, ...
        self_results(ip).batch_s,self_results(ip).speedup, ...
        self_results(ip).rel_inf);
end

function lambda = apply_one_body_loop(tau,U,Y,P,N_large,N_c)
PM = P*N_large;
lam_c_x = zeros(P*N_c,1);
lam_c_y = zeros(P*N_c,1);
U1 = U{1};
Y1 = Y{1};

for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c;
    tau_particle_x = tau((i-1)*N_large+1:i*N_large);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+i*N_large);

    lambda_i = Y1*(U1*[tau_particle_x; tau_particle_y]);
    lam_c_x(coarse_ind) = lambda_i(1:N_c);
    lam_c_y(coarse_ind) = lambda_i(N_c+1:2*N_c);
end

lambda = [lam_c_x; lam_c_y];
end

function lambda = apply_one_body_batched(tau,U,Y,P,N_large,N_c)
    PM = P*N_large;
    tau_x = reshape(tau(1:PM),N_large,P);
    tau_y = reshape(tau(PM+1:2*PM),N_large,P);
    tau_blocks = [tau_x; tau_y];
    
    step_blocks = U{1}*tau_blocks;
    lambda_blocks = Y{1}*step_blocks;
    
    lam_c_x = reshape(lambda_blocks(1:N_c,:),[],1);
    lam_c_y = reshape(lambda_blocks(N_c+1:2*N_c,:),[],1);
    lambda = [lam_c_x; lam_c_y];
end

function res = apply_self_block_loop(res,lambda_self,Nii,P,N_large,N_c)
    PM = P*N_large;
    n_coarse = P*N_c;
    lam_self_x = lambda_self(1:n_coarse);
    lam_self_y = lambda_self(n_coarse+1:2*n_coarse);
    
    for i = 1:P
        coarse_ind = (i-1)*N_c+1:i*N_c;
        tau_xy = [lam_self_x(coarse_ind); lam_self_y(coarse_ind)];
        uii = Nii*tau_xy;
    
        res((i-1)*N_large+1:i*N_large) = ...
            res((i-1)*N_large+1:i*N_large) - uii(1:end/2);
        res((i-1)*N_large+PM+1:i*N_large+PM) = ...
            res((i-1)*N_large+PM+1:i*N_large+PM) - uii(end/2+1:end);
    end
end

function res = apply_self_block_batched(res,lambda_self,Nii,P,N_large,N_c)
    PM = P*N_large;
    n_coarse = P*N_c;
    lambda_self_blocks = [
        reshape(lambda_self(1:n_coarse),N_c,P);
        reshape(lambda_self(n_coarse+1:2*n_coarse),N_c,P)
    ];
    
    uii_blocks = Nii*lambda_self_blocks;
    
    res(1:PM) = res(1:PM) - reshape(uii_blocks(1:N_large,:),[],1);
    res(PM+1:2*PM) = res(PM+1:2*PM) - ...
        reshape(uii_blocks(N_large+1:2*N_large,:),[],1);
end
