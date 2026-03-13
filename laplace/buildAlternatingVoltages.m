function v_body = buildAlternatingVoltages(q,R)
%BUILDALTERNATINGVOLTAGES Assign high-contrast body voltages from geometry.
%
% Syntax:
%   v_body = buildAlternatingVoltages(q,R)
%
% Inputs:
%   q - Complex particle centers (P x 1).
%   R - Scalar radius or per-body radii (P x 1).
%
% Output:
%   v_body - Per-body voltage values chosen to separate nearby bodies.
%
% Notes:
%   Uses a small set of alternating voltage levels and a greedy graph-based
%   assignment so near-neighbour discs tend to receive contrasting values.
%
% Anna Broms, Mar 2026

q = q(:);
P = numel(q);

if P==0
    v_body = zeros(0,1);
    return
end

if P==1
    v_body = 1;
    return
end

if isscalar(R)
    rad = repmat(R,P,1);
else
    rad = R(:);
end

D = abs(q-q.');
D(1:P+1:end) = inf;
gap_scale = rad + rad.';
score = D./gap_scale;
nearest = min(score,[],2);
thresh = 1.1*min(nearest);
adj = score <= thresh;
adj(1:P+1:end) = false;

levels = [1.5; -1.5; 0.5; -0.5];
v_body = nan(P,1);

deg = sum(adj,2);
[~,order] = sortrows([-deg, real(q), imag(q)]);

for idx = 1:P
    k = order(idx);
    nbr = find(adj(k,:) & ~isnan(v_body.'));

    if isempty(nbr)
        used = v_body(~isnan(v_body));
        if isempty(used)
            v_body(k) = levels(1);
        else
            score_level = abs(levels - mean(used));
            [~,best] = max(score_level);
            v_body(k) = levels(best);
        end
        continue
    end

    nbr_vals = v_body(nbr);
    min_sep = zeros(numel(levels),1);
    sum_sep = zeros(numel(levels),1);
    for j = 1:numel(levels)
        sep = abs(levels(j)-nbr_vals);
        min_sep(j) = min(sep);
        sum_sep(j) = sum(sep);
    end

    score_level = [min_sep, sum_sep, abs(levels)];
    best = pickBestLevel(score_level);
    v_body(k) = levels(best);
end

for sweep = 1:2
    for idx = 1:P
        k = order(idx);
        nbr = find(adj(k,:));
        if isempty(nbr)
            continue
        end

        nbr_vals = v_body(nbr);
        min_sep = zeros(numel(levels),1);
        sum_sep = zeros(numel(levels),1);
        for j = 1:numel(levels)
            sep = abs(levels(j)-nbr_vals);
            min_sep(j) = min(sep);
            sum_sep(j) = sum(sep);
        end

        score_level = [min_sep, sum_sep, abs(levels)];
        best = pickBestLevel(score_level);
        v_body(k) = levels(best);
    end
end

end

function best = pickBestLevel(score_level)
[~,ord] = sortrows([-score_level(:,1),-score_level(:,2),-score_level(:,3)]);
best = ord(1);
end
