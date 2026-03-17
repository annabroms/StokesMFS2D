function beta_rot = rotateStokesPairSourceVector(beta,Nf,im_i,im_j,phase_f,vec_rot)
%ROTATESTOKESPAIRSOURCEVECTOR Rotate pair-ordered fine/image Stokes sources.
%
% Ordering:
%   [f_i_x; e_i_x; f_j_x; e_j_x; f_i_y; e_i_y; f_j_y; e_j_y]

if nargin < 6 || isempty(vec_rot)
    vec_rot = 1;
end

beta = beta(:);
if isempty(beta)
    beta_rot = beta;
    return
end

f_ind1_x = 1:Nf;
e_ind1_x = Nf+1:Nf+im_i;
f_ind2_x = Nf+im_i+1:2*Nf+im_i;
e_ind2_x = 2*Nf+im_i+1:2*Nf+im_i+im_j;
f_ind1_y = 2*Nf+im_i+im_j+1:3*Nf+im_i+im_j;
e_ind1_y = 3*Nf+im_i+im_j+1:3*Nf+2*im_i+im_j;
f_ind2_y = 3*Nf+2*im_i+im_j+1:4*Nf+2*im_i+im_j;
e_ind2_y = 4*Nf+2*im_i+im_j+1:4*Nf+2*im_i+2*im_j;

fine_pair = [beta(f_ind1_x) beta(f_ind2_x) beta(f_ind1_y) beta(f_ind2_y)];
fine_pair = rotatePairOrderedStokesData(fine_pair,Nf,phase_f,vec_rot);

ei = vec_rot*(beta(e_ind1_x) + 1i*beta(e_ind1_y));
ej = vec_rot*(beta(e_ind2_x) + 1i*beta(e_ind2_y));

beta_rot = zeros(size(beta));
beta_rot(f_ind1_x) = fine_pair(:,1);
beta_rot(f_ind2_x) = fine_pair(:,2);
beta_rot(f_ind1_y) = fine_pair(:,3);
beta_rot(f_ind2_y) = fine_pair(:,4);
beta_rot(e_ind1_x) = real(ei);
beta_rot(e_ind1_y) = imag(ei);
beta_rot(e_ind2_x) = real(ej);
beta_rot(e_ind2_y) = imag(ej);
end
