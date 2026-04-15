function triplet = getStokesTripletInstance(triplet_cache,row)
%GETSTOKESTRIPLETINSTANCE Recover one actual Stokes triplet from a cached group.

meta = triplet_cache.meta(row);
group = triplet_cache.groups(meta.group_id);

triplet = struct();
triplet.meta = meta;
triplet.group = group;
triplet.q_trip = meta.anchor + meta.rot*group.q_trip;
triplet.rimage_cells = cell(3,1);
triplet.refine_cells = cell(3,1);
for body = 1:3
    triplet.rimage_cells{body} = map_points_from_canonical(group.rimage_canon{body},meta);
    triplet.refine_cells{body} = map_points_from_canonical(group.refine_canon{body},meta);
end

end

function z_actual = map_points_from_canonical(z_canon,meta)
z_canon = z_canon(:);
if isempty(z_canon)
    z_actual = zeros(0,1);
    return
end

z_actual = meta.anchor + meta.rot*z_canon;
end
