clear; close all

%What I learnt: With a goal of an independent number of iterations with P, 
% the number of singular value modes depends on delta. For the smallest
% delta tested, 0.1, k = 20 modes are sufficient. 


deltavec = [0.1, 0.5, 1, 5];

Pvec = [10 50 100 200];

svec = [1:10 20:10:50];
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
                [FT,lambda, it, gmres_tol, err] = solve_res_1B(q,U,W,rads,images,lr,visualise);
                it_lm(l,m) = it;
                err_lm(l,m) = err;
            end

        end
        iters{k,i} = it_lm;
        errs{k,i} = err_lm;
    end
end


%%
%save('data/oct18_deflation_Pdelta.mat','Pvec','errs',"iters","svec","deltavec");

c = {'b','r','k','c'};
m = {'o','s','d','+'};

for k = 1:length(Pvec)
    for i = 1:length(deltavec)
        it_lm = iters{k,i};
        it_plot(i,:) = it_lm;
        err = errs{k,i};

        figure(3);
        semilogy(svec,it_lm,'-','Color',c{k},'Marker',m{i});
        hold on

        figure(4);
        semilogy(svec,err,'-','Color',c{k},'Marker',m{i});
        hold on

        
    end
end
figure(3)
axis tight
grid on
xlabel('$k$','Interpreter','latex')
ylabel('GMRES iterations','Interpreter','latex')

%%
c = winter(length(svec));  
for i = 1:length(deltavec)
    for k = 1:length(Pvec)
        it_plot(k,:) = iters{k,i};
        


    end
    figure(5);
    subplot(1,4,i)
    for k = 1:length(svec)
        plot(Pvec,it_plot(:,k),'.-','Color',c(k,:));
        hold on
    end
    ylim([0,150])
    xlabel('$P$','Interpreter','latex','FontSize',16)
    ylabel('Iterations','Interpreter','latex','FontSize',16)
    title(['$\delta$ = ', num2str(deltavec(i))],'Interpreter','latex','FontSize',16)
    xlim([10,200])
    grid on
end

%%

