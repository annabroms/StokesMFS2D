function Ntrip = evaluateCoarseOnTriplet(q_trip,rbase_in_c,rout_cells)
%EVALUATECOARSEONTRIPLET Evaluate coarse cross-interactions on a triplet.
%
% Syntax:
%   Ntrip = evaluateCoarseOnTriplet(q_trip,rbase_in_c,rout_cells)
%
% Inputs:
%   q_trip      - 3x1 complex triplet centers
%   rbase_in_c  - coarse proxy locations relative to each center
%   rout_cells  - 3x1 cell array of target locations, one cell per body
%
% Output ordering:
%   sources: [x_1; x_2; x_3; y_1; y_2; y_3]
%   values : [u_x1; u_x2; u_x3; u_y1; u_y2; u_y3]

if nargin < 3
    error('evaluateCoarseOnTriplet requires q_trip, rbase_in_c, and rout_cells.');
end

q_trip = q_trip(:);
if numel(q_trip) ~= 3 || numel(rout_cells) ~= 3
    error('evaluateCoarseOnTriplet expects exactly three bodies.');
end

n_src = numel(rbase_in_c);
n_tgt = cellfun(@numel,rout_cells(:));
ntot = sum(n_tgt);

Nxx = zeros(ntot,3*n_src);
Nxy = zeros(ntot,3*n_src);
Nyx = zeros(ntot,3*n_src);
Nyy = zeros(ntot,3*n_src);

row_start = [1; 1 + n_tgt(1); 1 + n_tgt(1) + n_tgt(2)];
col_start = [1, 1 + n_src, 1 + 2*n_src];

for tgt_body = 1:3
    rt = rout_cells{tgt_body};
    row_idx = row_start(tgt_body):row_start(tgt_body) + n_tgt(tgt_body) - 1;
    for src_body = 1:3
        if src_body == tgt_body
            continue
        end

        Nbt = stokSLPmat(rbase_in_c + q_trip(src_body),rt,1);
        xrow = 1:n_tgt(tgt_body);
        yrow = n_tgt(tgt_body) + (1:n_tgt(tgt_body));
        xx = Nbt(xrow,1:n_src);
        xy = Nbt(xrow,n_src+1:end);
        yx = Nbt(yrow,1:n_src);
        yy = Nbt(yrow,n_src+1:end);

        col_idx = col_start(src_body):col_start(src_body) + n_src - 1;
        Nxx(row_idx,col_idx) = xx;
        Nxy(row_idx,col_idx) = xy;
        Nyx(row_idx,col_idx) = yx;
        Nyy(row_idx,col_idx) = yy;
    end
end

Ntrip = [Nxx, Nxy; Nyx, Nyy];

end
