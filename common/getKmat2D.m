function K = getKmat2D(r,q)
%GETKMAT2D(r,q) returns the matrix K that maps RBM to velocity at the boundary 
% for particles in 2D. The transpose, K', maps given force denisty to force and
% torque. Particles centered at q and discretised in points specified by the 
% (complex valued) positions in r.

cross_mat = @(x) [-imag(x); real(x)];

%for reference regarding the sign, see https://en.wikipedia.org/wiki/Cross_product 

NP = size(r,1); 
P = size(q,1); 
N = NP/P; %number of points per particle


%debug = 1;


r = transpose(r); 
q = transpose(q); 

K = zeros(2*NP,3*P);

%crete positive identity blocks
%J = repmat(eye(2),na,1);

%order x, then instead
J = blkdiag(ones(N,1),ones(N,1));


for i = 1:P %loop over particles
    B = zeros(2*N,1);
    %create particle block
    for k = 1:N
        %B(2*(k-1)+1:2*k,:) = cross_mat(r((i-1)*na+k)-q(i));  
        %B(2*(k-1)+1:2*k) = cross_mat(r((i-1)*na+k)-q(i));

        %very inefficient.
        C = cross_mat(r((i-1)*N+k)-q(i)); 
        B(k) = C(1);
        B(k+N) = C(2); 
    end
   
    K_block = [J,B];
    K((i-1)*2*N+1:i*2*N,3*(i-1)+1:3*i) = K_block;
end

end