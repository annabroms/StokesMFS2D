function K = getKmat2D(r,q)
%GETKMAT2D(r,q) returns the matrix K that maps a vector of RBM to velocities at the boundary 
% for particles in 2D. The transpose, K', maps a given force density to net force and
% torque. Particles centered at q and discretised in points specified by the 
% (complex valued) positions in r.

cross_mat = @(x) [-imag(x); real(x)]; %2D cross product

%for reference regarding the sign, see https://en.wikipedia.org/wiki/Cross_product 

NP = size(r,1); 
P = size(q,1); 
N = NP/P; %number of points per particle

r = r.';  
q = q.'; 

K = zeros(2*NP,3*P);

%crete "identity" blocks, sorted with first x, then y
J = blkdiag(ones(N,1),ones(N,1));

for i = 1:P %loop over particles
    B = zeros(2*N,1);
    %create particle block
    for k = 1:N %loop over sources on this particle

        %very inefficient...
        C = cross_mat(r((i-1)*N+k)-q(i)); 
        B(k) = C(1);
        B(k+N) = C(2); 
    end
   
    K_block = [J,B]; %first translational degrees of freedom, then rotational
    K((i-1)*2*N+1:i*2*N,3*(i-1)+1:3*i) = K_block;
end

end