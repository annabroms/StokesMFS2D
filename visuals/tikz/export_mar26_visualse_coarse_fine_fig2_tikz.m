clear;
close all;
clc;

% Convenience entry point: this regenerates the figure-2 PDF/demo files by
% calling the shared exporter, which currently writes both figure 2 and
% figure 3 assets in one pass.
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir,'export_mar26_visualse_coarse_fine_fig3_tikz.m'));
