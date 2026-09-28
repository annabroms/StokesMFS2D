%% SEP25_CAPACITANCE_CMAP_GAP_INTERP_SUBDOMAINS_ALPHA
% Repeat the equally subdivided Sep. 25 capacitance-Cmap interpolation
% experiment in alpha=acosh(1+delta/(2R)) instead of log(delta/R).
%
% All numerical settings, including Tikhonov regularisation, the constant
% 150-node ellipse construction, candidate local orders, validation grids,
% visible figures, and the two base gap ranges, are owned by the shared
% experiment script.  This wrapper changes only the interpolation
% coordinate and writes alpha-tagged result and figure files under data/.
%
% Anna Broms, Sep 25, 2026

cap_interp_coordinate_override = 'alpha';
run(fullfile(fileparts(mfilename('fullpath')), ...
    'sep25_capacitance_cmap_gap_interp_subdomains.m'));
