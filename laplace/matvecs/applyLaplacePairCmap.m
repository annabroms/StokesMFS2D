function [tau_pair,pair_qv] = applyLaplacePairCmap( ...
        pair_cache,group,rhs_pair)
%APPLYLAPLACEPAIRCMAP Apply an exact, full, or factorised interpolated map.

map_kind = getOptField(group,'map_kind','none');
switch map_kind
    case {'','none','full'}
        tau_pair = group.Cmap*rhs_pair;
    case 'reduced'
        panel = pair_cache.interpolator.panels(group.interp_panel);
        tau_pair = panel.Cref*rhs_pair + ...
            panel.U*(group.interp_B*(panel.V'*rhs_pair));
    case 'reduced_noconst'
        panel = pair_cache.interpolator.panels(group.interp_panel);
        tau_pair = panel.U*(group.interp_B*(panel.V'*rhs_pair));
    otherwise
        error('applyLaplacePairCmap:BadMapKind', ...
            'Unknown pair map kind "%s".',map_kind);
end
pair_qv = group.Cmap_QV*rhs_pair;
end
