function panels = buildLaplaceCmapAlphaPanelPlan(opt)
%BUILDLAPLACECMAPALPHAPANELPLAN Reproducible calibrated alpha-panel layout.
%
% Four equal-alpha panels cover [delta_min,delta_star], all with q=5.
% Eight equal-alpha panels cover [delta_star,delta_pair].  The documented
% q=7 fifth panel is retained, and the fourth is also assigned q=7 because
% the production interval includes delta_pair rather than 0.999*delta_pair.
% Cmap_QV is only 2-by-160, so it uses an independent q=17 order on the
% first small-gap panel and q=9 elsewhere.  This removes the charge-map
% error without lengthening the expensive Cmap or B interpolation sums.

R = opt.rad;
delta_min = opt.smallest_delta;
delta_max = opt.delta_pair;
delta_star = (R-opt.Rp_f)^2/opt.Rp_f;
if ~(delta_min < delta_star && delta_star < delta_max)
    error('buildLaplaceCmapAlphaPanelPlan:BadRanges', ...
        'Expected smallest_delta < delta_star < delta_pair.');
end

counts = [4 8];
q_by_base = {[5 5 5 5],[5 5 5 7 7 5 5 5]};
q_charge_by_base = {[17 9 9 9],9*ones(1,8)};
rank_by_base = {[44 46 49 48],[52 47 54 45 42 43 44 44]};
range_delta = [delta_min delta_star; delta_star delta_max];
to_alpha = @(delta) acosh(1+delta/(2*R));
to_delta = @(alpha) 2*R*(cosh(alpha)-1);

template = struct('base_index',[],'local_index',[],'alpha_lo',[], ...
    'alpha_hi',[],'delta_lo',[],'delta_hi',[],'q',[],'q_charge',[], ...
    'spatial_rank',[],'alpha_nodes',[],'delta_nodes',[], ...
    'barycentric_weights',[],'reference_delta',[], ...
    'charge_alpha_nodes',[],'charge_delta_nodes',[], ...
    'charge_barycentric_weights',[],'basis_alpha_nodes',[], ...
    'basis_delta_nodes',[]);
panels = repmat(template,sum(counts),1);
index = 0;
for ib = 1:2
    alpha_range = to_alpha(range_delta(ib,:));
    edges = linspace(alpha_range(1),alpha_range(2),counts(ib)+1);
    for local_index = 1:counts(ib)
        index = index+1;
        q = q_by_base{ib}(local_index);
        q_charge = q_charge_by_base{ib}(local_index);
        alpha_lo = edges(local_index);
        alpha_hi = edges(local_index+1);
        alpha_nodes = chebLobatto(alpha_lo,alpha_hi,q);
        panels(index).base_index = ib;
        panels(index).local_index = local_index;
        panels(index).alpha_lo = alpha_lo;
        panels(index).alpha_hi = alpha_hi;
        panels(index).delta_lo = to_delta(alpha_lo);
        panels(index).delta_hi = to_delta(alpha_hi);
        panels(index).q = q;
        panels(index).q_charge = q_charge;
        panels(index).spatial_rank = rank_by_base{ib}(local_index);
        panels(index).alpha_nodes = alpha_nodes;
        panels(index).delta_nodes = to_delta(alpha_nodes);
        panels(index).barycentric_weights = barycentricWeights(q);
        panels(index).charge_alpha_nodes = ...
            chebLobatto(alpha_lo,alpha_hi,q_charge);
        panels(index).charge_delta_nodes = ...
            to_delta(panels(index).charge_alpha_nodes);
        panels(index).charge_barycentric_weights = ...
            barycentricWeights(q_charge);
        panels(index).reference_delta = ...
            (panels(index).delta_lo+panels(index).delta_hi)/2;
        panels(index).basis_alpha_nodes = ...
            chebLobatto(alpha_lo,alpha_hi,9);
        panels(index).basis_delta_nodes = ...
            to_delta(panels(index).basis_alpha_nodes);
    end
end
end

function nodes = chebLobatto(lo,hi,q)
j = (0:q-1)';
nodes = (lo+hi)/2+(hi-lo)/2*cos(j*pi/(q-1));
end

function weights = barycentricWeights(q)
j = (0:q-1)';
weights = (-1).^j;
weights([1 end]) = weights([1 end])/2;
end
