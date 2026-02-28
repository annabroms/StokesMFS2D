%test lr preconditioning

%What I learnt: Smaller delta requires a larger number of modes for a small
%number of iterations. The same test script was used both to test the
%performance with trigonometric modes and with modes corresponding to the
%first k singular vectors of the single particle SVD. 

close all; clear;
images = 0; %images not needed for well separated particles
deltavec = [0.1:0.1:1 2 5]; 
rng(8); 

P = 10;
nruns = 10; 
rads = ones(P,1);
visualise = 0; 
lrvec = [0 4:10 23:10:70];
lrvec = [0 5 7 11:4:27 (20:10:30)*4+3]; 
% SVD k chosen to match the number of coarse basis functions used with a trigonometric basis. should be up to 60 here

for k = 1:length(deltavec)
    delta = deltavec(k);
    for i = 1:nruns
        [q,B] = grow_cluster(P,delta,2);
        U = rand(P,2); %translational velocities 
        W = rand(P,1); %angular velocities 
        for l = 1:length(lrvec)
            lr = lrvec(l);
            [FT,lambda, it, gmres_tol, err] = solve_res_1B(q,U,W,rads,images,lr,visualise);
            if it > 200
                disp('Waiting')
            end
            iters(k,i,l) = it;
            errs(k,i,l) = err;
        end
    end
end

%%
% Visualise the results
c = lines(length(lrvec)+1);
figure(4)
loglog(nan,nan);
hold on

figure(3)
semilogx(nan,nan)
hold on

markers = {'*','.','s','d','o'};

for i = 1:length(lrvec)
    
    figure(1)
    %loglog(nan,nan);
    %hold on
    upper_bound = max(reshape(iters(:,:,i),length(deltavec),nruns)');
    lower_bound = min(reshape(iters(:,:,i),length(deltavec),nruns)');
    fill([deltavec, fliplr(deltavec)], [upper_bound, fliplr(lower_bound)], c(i,:), 'FaceAlpha', 0.3);
    hold on
    
    figure(2)
    semilogy(nan,nan);
    hold on
    upper_bound = max(reshape(errs(:,:,i),length(deltavec),nruns)');
    lower_bound = min(reshape(errs(:,:,i),length(deltavec),nruns)');
    fill([deltavec, fliplr(deltavec)], [upper_bound, fliplr(lower_bound)],c(i,:), 'FaceAlpha', 0.3);
    hold on

    
    figure(3)
    mean_val = mean(reshape(iters(:,:,i),length(deltavec),nruns)');
    if i>8
        h3(i) = plot(deltavec,mean_val,'-','Marker',markers{mod(i+1,5)+1},'Color',c(i+1,:));
    else
        h3(i) = plot(deltavec,mean_val,'-','Marker',markers{mod(i,5)+1},'Color',c(i,:));
    end
    

    figure(4)
    mean_val = mean(reshape(errs(:,:,i),length(deltavec),nruns)');
    if i>8
        h4(i) = plot(deltavec,mean_val,'-','Marker',markers{mod(i+1,5)+1},'Color',c(i+1,:));
    else
        h4(i) = plot(deltavec,mean_val,'-','Marker',markers{mod(i,5)+1},'Color',c(i,:));
    end
    
    

end

figure(3)
axis tight
xlabel('$\delta$','interpreter','latex')
ylabel('GMRES iterations','interpreter','latex')
grid on

labels = arrayfun(@(l) sprintf('$k = %d$', l), [lrvec(2:end)-3], 'UniformOutput', false);

%labels = arrayfun(@(l) sprintf('$\\text{mode } %d$', l), [lrvec(2:end)-3], 'UniformOutput', false);
legend(h3,['No deflation', labels], 'Interpreter', 'latex', 'Location', 'best');

figure(4)
axis tight
xlabel('$\delta$','interpreter','latex')
ylabel('Relative surface residual','interpreter','latex')
grid on
legend(h4,['No deflation', labels], 'Interpreter', 'latex', 'Location', 'best');

%save('data/oct15_deflation_sweep_singvals.mat','iters','errs','deltavec','P','lrvec','nruns')