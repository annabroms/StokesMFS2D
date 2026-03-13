clear;
close all;
set(0,'DefaultFigureVisible','off');

startup;
rng(13);

ns = 60;
nt = 72;
nrepeat = 30000;

rsrc = rand(ns,1) + 1i*rand(ns,1);
rtar = rand(nt,1) + 1i*rand(nt,1);
sigma = randn(ns,1);

fprintf('=== Laplace Direct vs r2ddir Timing Experiment (Mar 13, 2026) ===\n');
fprintf('sources=%d, targets=%d, repeats=%d\n',ns,nt,nrepeat);

if exist('r2ddir','file')~=2
    error('mar13_laplace_direct_vs_r2ddir:r2ddirMissing', ...
        'r2ddir is not on path. Run startup.m and make sure ../fmm2d is available.');
end

u_matlab = lapSLPdirect(rsrc,rtar,sigma);
u_r2ddir = eval_lap_slp_r2ddir(rsrc,rtar,sigma);

abs_err = norm(u_r2ddir-u_matlab,inf);
rel_err = abs_err/max(1,norm(u_matlab,inf));

fprintf('warm-up check: abs err %.3e, rel err %.3e\n',abs_err,rel_err);

tic;
for k = 1:nrepeat
    u_matlab = lapSLPdirect(rsrc,rtar,sigma);
end
t_matlab = toc;

tic;
for k = 1:nrepeat
    u_r2ddir = eval_lap_slp_r2ddir(rsrc,rtar,sigma);
end
t_r2ddir = toc;

fprintf('lapSLPdirect total time: %.6f s\n',t_matlab);
fprintf('r2ddir      total time: %.6f s\n',t_r2ddir);
fprintf('lapSLPdirect time/call: %.3e s\n',t_matlab/nrepeat);
fprintf('r2ddir      time/call: %.3e s\n',t_r2ddir/nrepeat);
fprintf('speedup (lapSLPdirect / r2ddir): %.2fx\n',t_matlab/max(t_r2ddir,eps));

function u = eval_lap_slp_r2ddir(rsrc,rtar,sigma)
srcinfo = struct();
srcinfo.sources = [real(rsrc)'; imag(rsrc)'];
srcinfo.nd = size(sigma,2);
srcinfo.charges = (-sigma.')/(2*pi);

targ = [real(rtar)'; imag(rtar)'];
U = r2ddir(srcinfo,targ,1);
u = U.pottarg.';
end
