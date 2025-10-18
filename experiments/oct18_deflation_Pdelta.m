clear; close all

deltavec = [0.1, 0.5, 1, 5];

Pvec = [10 50 100 200];

svec = [1:10 20:10:100];
nruns = 1;
rng(8);
visualise = 0; 
images = 0; 

for k = 1:length(Pvec)
    P = Pvec(k);
    rads = ones(P,1);
    for i = 1:length(deltavec)
        delta = deltavec(i);

        for l = 1:nruns
            %generate geometry here
           
            [q,B] = grow_cluster(P,delta,2);
            
            U = rand(P,2); %translational velocities 
            W = rand(P,1); %angular velocities 
            
            for m = 1:length(svec)
                lr = svec(m)+3;
                [FT,lambda, it, gmres_tol, err] = solve_2D_res(q,U,W,rads,images,lr,visualise);
                it_lm(l,m) = it;
                err_lm(l,m) = err;
            end

        end
        iters{k,i} = it_lm;
        errs{k,i} = err;
    end
end
