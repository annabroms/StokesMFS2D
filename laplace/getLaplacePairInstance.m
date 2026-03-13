function pair = getLaplacePairInstance(pair_cache,row)
%GETLAPLACEPAIRINSTANCE Recover one actual pair from a canonical cached group.

meta = pair_cache.meta(row);
group = pair_cache.groups(meta.group_id);

pair = struct();
pair.meta = meta;
pair.group = group;
pair.rimage_i = map_points_from_canonical(group.rimage_canon{1},meta);
pair.rimage_j = map_points_from_canonical(group.rimage_canon{2},meta);

end

function z_actual = map_points_from_canonical(z_canon,meta)
z_canon = z_canon(:);
if isempty(z_canon)
    z_actual = zeros(0,1);
    return
end

z_actual = meta.mid + meta.rot*z_canon;

end
