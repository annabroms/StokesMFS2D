function beta_rot = rotateStokesTripletSourceVector(beta,Nf,im_counts,phase_f,vec_rot)
%ROTATESTOKESTRIPLETSOURCEVECTOR Rotate triplet fine/image Stokes sources.
%
% Ordering:
%   [f1x; e1x; f2x; e2x; f3x; e3x; f1y; e1y; f2y; e2y; f3y; e3y]

if nargin < 5 || isempty(vec_rot)
    vec_rot = 1;
end

beta = beta(:);
if isempty(beta)
    beta_rot = beta;
    return
end

idx = get_triplet_indices(Nf,im_counts);
fine_trip = [beta(idx.f1x), beta(idx.f2x), beta(idx.f3x), ...
             beta(idx.f1y), beta(idx.f2y), beta(idx.f3y)];
fine_trip = rotatePairOrderedStokesData(fine_trip,Nf,phase_f,vec_rot);

e1 = vec_rot*(beta(idx.e1x) + 1i*beta(idx.e1y));
e2 = vec_rot*(beta(idx.e2x) + 1i*beta(idx.e2y));
e3 = vec_rot*(beta(idx.e3x) + 1i*beta(idx.e3y));

beta_rot = zeros(size(beta));
beta_rot(idx.f1x) = fine_trip(:,1);
beta_rot(idx.f2x) = fine_trip(:,2);
beta_rot(idx.f3x) = fine_trip(:,3);
beta_rot(idx.f1y) = fine_trip(:,4);
beta_rot(idx.f2y) = fine_trip(:,5);
beta_rot(idx.f3y) = fine_trip(:,6);
beta_rot(idx.e1x) = real(e1);
beta_rot(idx.e1y) = imag(e1);
beta_rot(idx.e2x) = real(e2);
beta_rot(idx.e2y) = imag(e2);
beta_rot(idx.e3x) = real(e3);
beta_rot(idx.e3y) = imag(e3);

end

function idx = get_triplet_indices(Nf,im_counts)
im1 = im_counts(1);
im2 = im_counts(2);
im3 = im_counts(3);

idx = struct();
idx.f1x = 1:Nf;
idx.e1x = Nf+1:Nf+im1;
idx.f2x = Nf+im1+1:2*Nf+im1;
idx.e2x = 2*Nf+im1+1:2*Nf+im1+im2;
idx.f3x = 2*Nf+im1+im2+1:3*Nf+im1+im2;
idx.e3x = 3*Nf+im1+im2+1:3*Nf+im1+im2+im3;

base_y = 3*Nf + im1 + im2 + im3;
idx.f1y = base_y + (1:Nf);
idx.e1y = base_y + Nf + (1:im1);
idx.f2y = base_y + Nf + im1 + (1:Nf);
idx.e2y = base_y + 2*Nf + im1 + (1:im2);
idx.f3y = base_y + 2*Nf + im1 + im2 + (1:Nf);
idx.e3y = base_y + 3*Nf + im1 + im2 + (1:im3);
idx.ntot = 6*Nf + 2*(im1 + im2 + im3);
end
