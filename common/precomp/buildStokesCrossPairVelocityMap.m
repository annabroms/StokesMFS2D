function Ucross = buildStokesCrossPairVelocityMap(Epair,N_src,N_tgt)
%BUILDSTOKESCROSSPAIRVELOCITYMAP Keep only cross-body Stokeslet blocks.

Ucross = zeros(size(Epair));

tgt_i_x = 1:N_tgt;
tgt_j_x = N_tgt+1:2*N_tgt;
tgt_i_y = 2*N_tgt+1:3*N_tgt;
tgt_j_y = 3*N_tgt+1:4*N_tgt;

src_i_x = 1:N_src;
src_j_x = N_src+1:2*N_src;
src_i_y = 2*N_src+1:3*N_src;
src_j_y = 3*N_src+1:4*N_src;

Ucross(tgt_i_x,[src_j_x src_j_y]) = Epair(tgt_i_x,[src_j_x src_j_y]);
Ucross(tgt_j_x,[src_i_x src_i_y]) = Epair(tgt_j_x,[src_i_x src_i_y]);
Ucross(tgt_i_y,[src_j_x src_j_y]) = Epair(tgt_i_y,[src_j_x src_j_y]);
Ucross(tgt_j_y,[src_i_x src_i_y]) = Epair(tgt_j_y,[src_i_x src_i_y]);

Ucross = -Ucross;
end
