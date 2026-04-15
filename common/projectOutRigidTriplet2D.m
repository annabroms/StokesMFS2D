function triplet_proj = projectOutRigidTriplet2D(triplet_vec, rsrc_cells, q_trip)
%PROJECTOUTRIGIDTRIPLET2D Project each body of a triplet independently.
%
% Input ordering:
%   [x_1; x_2; x_3; y_1; y_2; y_3]

if nargin < 3
    error('projectOutRigidTriplet2D requires triplet_vec, rsrc_cells, and q_trip.');
end

if numel(rsrc_cells) ~= 3 || numel(q_trip) ~= 3
    error('projectOutRigidTriplet2D expects three source sets and three centers.');
end

triplet_vec = triplet_vec(:);
q_trip = q_trip(:);

n = cellfun(@numel,rsrc_cells(:));
if numel(triplet_vec) ~= 2*sum(n)
    error('projectOutRigidTriplet2D:SizeMismatch', ...
        'triplet_vec has the wrong length for the supplied source sets.');
end

x_blocks = split_blocks(triplet_vec(1:sum(n)),n);
y_blocks = split_blocks(triplet_vec(sum(n)+1:end),n);

triplet_proj = zeros(size(triplet_vec));
offset_x = 0;
offset_y = sum(n);
for body = 1:3
    tau_body = [x_blocks{body}; y_blocks{body}];
    tau_proj = projectOutRigid2D(tau_body,rsrc_cells{body},q_trip(body));
    nb = n(body);
    triplet_proj(offset_x + (1:nb)) = tau_proj(1:nb);
    triplet_proj(offset_y + (1:nb)) = tau_proj(nb+1:end);
    offset_x = offset_x + nb;
    offset_y = offset_y + nb;
end

end

function blocks = split_blocks(vec,counts)
blocks = cell(numel(counts),1);
offset = 0;
for k = 1:numel(counts)
    blocks{k} = vec(offset + (1:counts(k)));
    offset = offset + counts(k);
end
end
