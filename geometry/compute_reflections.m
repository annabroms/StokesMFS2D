function [z,zvec,itr] = compute_reflections(z0, r1, r2, tol,iter)


if nargin<4
    tol = 1e-15;
    iter = 0;
elseif nargin <5
    iter = 0; 
end

zvec = [z0];
t0 = z0;
 
if iter
    for k = 1:iter
        t1 = conj(r2(t0));
        t2 = conj(r1(t1));
        t0 = t2;
        zvec = [zvec; t0];
    end
    
else

    h = 1;
    itr = 0; 
    %tol = 1e-15; % How accurate are the AAA approximations of the curves?

    while h > tol
        t1 = conj(r2(t0));
        t2 = conj(r1(t1));
        h = abs(t0-t2); 
        t0 = t2;
        zvec = [zvec; t0];
        itr = itr+1; 
    end

end

z = zvec(end);



