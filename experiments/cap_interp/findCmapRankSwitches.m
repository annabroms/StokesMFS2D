function [switch_deltas,scan] = findCmapRankSwitches(delta_min,delta_max,R,opt,grids,n_scan,cache_file)
%FINDCMAPRANKSWITCHES Delta locations where the underlying pseudoinverse
% ranks or active-source counts of C(delta) change.
%
% [switch_deltas,scan] = findCmapRankSwitches(delta_min,delta_max,R,opt, ...
%   grids,n_scan,cache_file)
%
% Densely scans n_scan log-spaced gaps in [delta_min,delta_max] (the same
% diagnostic sep24_capacitance_cmap_delta_smoothness.m uses: rank_fine,
% rank_peanut from the two pseudoinverses inside getPairBasisLaplace, and
% n_image, the number of active ellipse-enhancement source nodes per
% body) to find candidate switch brackets, then bisects each bracket to a
% tight relative tolerance. This matters because Chebyshev-Lobatto nodes
% cluster arbitrarily close to subinterval endpoints as q grows: a switch
% location known only to the coarse scan's resolution (~n_scan points over
% the whole range) can still fall strictly inside a "clean" subinterval
% once refined, silently corrupting every training snapshot on that
% subinterval (global polynomial interpolation is not local: one wrong
% snapshot spoils the whole fit, not just nearby evaluations).
%
% switch_deltas collects, for every detected change, the delta value
% immediately after it (bisected to rel_tol below), so that
% [delta_min, switch_deltas(1)], [switch_deltas(1), switch_deltas(2)], ...,
% [switch_deltas(end), delta_max] are exactly the subintervals on which none
% of rank_fine/rank_peanut/n_image changes.
%
% scan is a struct with fields delta (n_scan x 1), rank_fine, rank_peanut,
% n_image (n_scan x 2), from the initial coarse scan, for further
% inspection/plotting.
%
% Anna Broms, Sep 2026

delta_scan = logspace(log10(delta_min),log10(delta_max),n_scan)';
[~,info_scan] = getOrBuildCmapSnapshots(delta_scan,R,opt,grids,cache_file);

scan.delta = delta_scan;
scan.rank_fine = [info_scan.rank_fine]';
scan.rank_peanut = [info_scan.rank_peanut]';
scan.n_image = reshape([info_scan.n_image],2,n_scan)';

sig = [scan.rank_fine,scan.rank_peanut,scan.n_image];
changed = find(any(diff(sig),2));

rel_tol = 1e-9;
switch_deltas = zeros(numel(changed),1);
for k = 1:numel(changed)
    i = changed(k);
    switch_deltas(k) = bisectRankSwitch(delta_scan(i),delta_scan(i+1),R,opt,grids,cache_file,rel_tol);
end

end

function delta_new = bisectRankSwitch(delta_old,delta_new,R,opt,grids,cache_file,rel_tol)
%BISECTRANKSWITCH Refine [delta_old,delta_new] (known to bracket exactly one
% signature change, delta_old on the old side) to the smallest delta_new
% (to relative precision rel_tol) at which the new signature is seen.
[~,info_old] = getOrBuildCmapSnapshots(delta_old,R,opt,grids,cache_file);
sig_old = signature(info_old);
max_iter = 60;
for it = 1:max_iter
    if (delta_new-delta_old) <= rel_tol*delta_new
        break
    end
    delta_mid = sqrt(delta_old*delta_new); % geometric mean: bisection in log(delta)
    [~,info_mid] = getOrBuildCmapSnapshots(delta_mid,R,opt,grids,cache_file);
    if isequal(signature(info_mid),sig_old)
        delta_old = delta_mid;
    else
        delta_new = delta_mid;
    end
end
end

function s = signature(info)
s = [info.rank_fine,info.rank_peanut,info.n_image];
end
