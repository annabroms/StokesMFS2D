%Test out the long range preconditioning

% Test of a constant basis in the x and y direction (lr = 1), or the rigid
% body matrix for each body (lr = 2). What I learnt: works well for dilute
% systems. Smaller reduction in the number of iterations for denser
% systems. 

clear;
close all;

Pvec = 10:10:100;
Pvec = [50 100];
Pvec = 200; 
%Pvec = Pvec(end-1); 
%Pvec = [10 40];
deltavec = [0.5, 1, 5, 10];
%deltavec = [5 10];
Nruns = 10;

visualise = 0;
images = 0;

for i = 1:length(deltavec)
    delta = deltavec(i);
    for k = 1:length(Pvec)
        P = Pvec(k);
        rad = ones(P,1);
        for j = 1:Nruns
            %generate a random configuration
            [q,B] = grow_cluster(P,delta,2);
            U = rand(P,2); %translational velocities 
            W = rand(P,1); %angular velocities 
            % lr = 0;
            % [~,~, it_0, ~, err_0] = solve_res_1B(q,U,W,rad,images,lr,visualise);
            % lr = 1;
            % [~,~, it_1, ~, err_1] = solve_res_1B(q,U,W,rad,images,lr,visualise);
            lr = 2;
            [~,~, it_2, ~, err_2] = solve_res_1B(q,U,W,rad,images,lr,visualise);

            iters_0(i,k,j) = it_0;
            iters_1(i,k,j) = it_1;
            iters_2(i,k,j) = it_2;

            errs_0(i,k,j) = err_0;
            errs_1(i,k,j) = err_1;
            errs_2(i,k,j) = err_2;
            
        end
    end
end
%%
for i = 1:length(deltavec)
    figure(1)
    plot(Pvec,reshape(errs_0(i,:,:),length(Pvec),Nruns),'*-')
    hold on
    plot(Pvec,reshape(errs_1(i,:,:),length(Pvec),Nruns),'o--')
    plot(Pvec,reshape(errs_2(i,:,:),length(Pvec),Nruns),'s:')

    figure(2)
    plot(Pvec,reshape(iters_0(i,:,:),length(Pvec),Nruns),'*-')
    hold on
    plot(Pvec,reshape(iters_1(i,:,:),length(Pvec),Nruns),'o--')
    plot(Pvec,reshape(iters_2(i,:,:),length(Pvec),Nruns),'s:')

end

%%

for i = 1:length(deltavec)
    figure(1)
    subplot(1,4,i)
    %loglog(nan,nan);
    %hold on
    upper_bound = max(reshape(iters_0(i,:,:),length(Pvec),Nruns)');
    lower_bound = min(reshape(iters_0(i,:,:),length(Pvec),Nruns)');
    fill([Pvec, fliplr(Pvec)], [upper_bound, fliplr(lower_bound)], 'r', 'FaceAlpha', 0.3);
    hold on
    upper_bound = max(reshape(iters_1(i,:,:),length(Pvec),Nruns)');
    lower_bound = min(reshape(iters_1(i,:,:),length(Pvec),Nruns)');
    fill([Pvec, fliplr(Pvec)], [upper_bound, fliplr(lower_bound)], 'b', 'FaceAlpha', 0.3);

    upper_bound = max(reshape(iters_2(i,:,:),length(Pvec),Nruns)');
    lower_bound = min(reshape(iters_2(i,:,:),length(Pvec),Nruns)');
    fill([Pvec, fliplr(Pvec)], [upper_bound, fliplr(lower_bound)], 'k', 'FaceAlpha', 0.3);
    %loglog(gapvec,1e-7*ones(size(gapvec)),'k--')
    ylim([0,100])
    xlabel('$P$','interpreter','latex')
    if i == 1
        ylabel('GMRES iterations','interpreter','latex')
    end
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on
    title(['$\delta=$ ', num2str(deltavec(i))], 'Interpreter', 'latex');
    xlim([min(Pvec),max(Pvec)])

    figure(2)
    subplot(1,4,i)
    %loglog(nan,nan);
    %hold on
    upper_bound = max(reshape(errs_0(i,:,:),length(Pvec),Nruns)');
    lower_bound = min(reshape(errs_0(i,:,:),length(Pvec),Nruns)');
    fill([Pvec, fliplr(Pvec)], [upper_bound, fliplr(lower_bound)], 'r', 'FaceAlpha', 0.3);
    hold on
    upper_bound = max(reshape(errs_1(i,:,:),length(Pvec),Nruns)');
    lower_bound = min(reshape(errs_1(i,:,:),length(Pvec),Nruns)');
    fill([Pvec, fliplr(Pvec)], [upper_bound, fliplr(lower_bound)], 'b', 'FaceAlpha', 0.3);

    upper_bound = max(reshape(errs_2(i,:,:),length(Pvec),Nruns)');
    lower_bound = min(reshape(errs_2(i,:,:),length(Pvec),Nruns)');
    fill([Pvec, fliplr(Pvec)], [upper_bound, fliplr(lower_bound)], 'k', 'FaceAlpha', 0.3);
    %loglog(gapvec,1e-7*ones(size(gapvec)),'k--')
   % ylim([0,100])
    xlabel('$P$','interpreter','latex')
    if i == 1
        ylabel('surface residual','interpreter','latex')
    end
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on
    title(['$\delta=$ ', num2str(deltavec(i))], 'Interpreter', 'latex');
    set(gca, 'YScale', 'log')
    ylim([1e-8, 1e-4])
    xlim([min(Pvec),max(Pvec)])
end