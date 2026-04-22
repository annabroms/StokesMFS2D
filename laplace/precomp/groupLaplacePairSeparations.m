function [group_id,group_sep,rep_rows] = groupLaplacePairSeparations(meta,sep_tol)
%GROUPLAPLACEPAIRSEPARATIONS Group close pairs by center separation.

sep = reshape([meta.sep],[],1);
[sep_sorted,order] = sort(sep);

group_id_sorted = zeros(size(sep_sorted));
group_sep = zeros(0,1);
rep_rows = zeros(0,1);

if isempty(sep_sorted)
    group_id = zeros(0,1);
    return
end

group_count = 1;
group_id_sorted(1) = group_count;
group_sep(1,1) = sep_sorted(1);
rep_rows(1,1) = order(1);

for k = 2:numel(sep_sorted)
    if abs(sep_sorted(k)-group_sep(group_count)) <= sep_tol
        group_id_sorted(k) = group_count;
    else
        group_count = group_count + 1;
        group_id_sorted(k) = group_count;
        group_sep(group_count,1) = sep_sorted(k);
        rep_rows(group_count,1) = order(k);
    end
end

group_id = zeros(size(sep));
group_id(order) = group_id_sorted;
end
