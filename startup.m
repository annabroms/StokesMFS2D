%STARTUP Add project paths and check external dependencies.
%
% Anna Broms, Mar 2026

repo_root = fileparts(mfilename('fullpath'));
addpath(genpath(repo_root));

% Expected location of memorygraph relative to this repository.
memorygraph_path = fullfile(repo_root,'..','memorygraph');
if isfolder(memorygraph_path)
    addpath(genpath(memorygraph_path));
    if exist('memorygraph','file') ~= 2
        fprintf(2,'[startup] memorygraph is not on path after startup.\n');
        fprintf(2,'[startup] opt.RAM_check requires memorygraph.\n');
    end
else
    fprintf(2,'[startup] memorygraph was not found at:\n  %s\n',memorygraph_path);
    fprintf(2,'[startup] Suggested fix:\n');
    fprintf(2,'  1) Clone https://github.com/ahbarnett/memorygraph\n');
    fprintf(2,'  2) Place it at the path above (or edit startup.m -> memorygraph_path)\n');
    fprintf(2,'[startup] opt.RAM_check requires memorygraph.\n');
end

% Expected location of flatiron fmm2d relative to this repository.
fmm2d_path = fullfile(repo_root,'..','fmm2d');
if isfolder(fmm2d_path)
    addpath(genpath(fmm2d_path));
else
    fprintf(2,'[startup] fmm2d was not found at:\n  %s\n',fmm2d_path);
    fprintf(2,'[startup] Suggested fix:\n');
    fprintf(2,'  1) Clone https://github.com/flatironinstitute/fmm2d\n');
    fprintf(2,'  2) Place it at the path above (or edit startup.m -> fmm2d_path)\n');
end

if exist('rfmm2d','file')~=3 && exist('rfmm2d','file')~=2
    fprintf(2,'[startup] rfmm2d is not on path after startup.\n');
    fprintf(2,'[startup] If fmm2d exists, compile/build its MATLAB interface and re-run startup.\n');
    fprintf(2,'[startup] Laplace/Stokes solves can fall back to direct summation with opt.use_fmm = 0, but will be much slower and is not recommended.\n');
end

