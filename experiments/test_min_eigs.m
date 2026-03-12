clear; 
close all; 

%test min eigs
deltavec = [1e-3 1e-2 1e-1];

N = 40;
%alpha_vec = linspace(0,2*pi/3,N);
alpha_vec = linspace(7*pi/12,2*pi/3,N);  
alpha_vec = linspace(2,2*pi/3,N); 
delta_pair = 5; %all pairs close to each other

delta_vec2 = logspace(-3,0,N);

t = linspace(0,2*pi);

fig1 = figure(1);
fig2 = figure(2);
fig3 = figure(3);
fig4 = figure(4);

test_name = 'triange_alphasweep_noimages_Nf1000.mat';




for l = 1:3
    %delta = 0.001;
    delta = deltavec(l);
    
    q = [0; 2+delta];
    
    for k = 1:N
        %load(test_name);
        alpha = alpha_vec(k); 
        %delta2 = delta_vec2(k); 
        q_all = [q; q(2)+(2+delta)*(cos(alpha)+1i*sin(alpha))];
        
    
        %q_all = [q; q(2)+(2+delta)*(cos(2*pi/3)+1i*sin(2*pi/3))+delta2*1i];
        dist(k,l) = abs(q_all(1)-q_all(3))-2;
    
    
    %     figure(3)
    %     for i = 1:3
    %         plot(real(q_all(i))+cos(t),imag(q_all(i))+sin(t))
    %         hold on
    %     end
    
        %keep track of only the relevant figures. 

        allFigures = findall(0, 'Type', 'figure');

        allFigures(allFigures == fig1) = []; % Remove fig1 from the list
        allFigures(allFigures == fig2) = []; 
        allFigures(allFigures == fig3) = []; 
        allFigures(allFigures == fig4) = []; 
        close(allFigures);

        [min_eig,c,noimage] = get_eigs(q_all,delta_pair);
        conds(k,l) = c;
        eigmins(k,l) = min_eig; 
        image_13(k,l) = noimage;

        save(test_name,'conds','alpha_vec','deltavec','eigmins','image_13');
    
    end
    


    figure(1)
    semilogy(alpha_vec,conds,'+-')
    %semilogy(delta_vec2,conds,'+-')
    hold on
    grid on
    axis tight
    %xticks([pi/2, 2*pi/3]);
    %xticks([7*pi/12, 2*pi/3]);
    %xticklabels({'7pi/12', '2\pi/3'});
    xlabel('$\alpha$','interpreter','latex')
   % xlabel('$\delta_2$','Interpreter','latex')
    ylabel('Condition number','Interpreter','latex')
    
    figure(2)
    semilogy(alpha_vec,eigmins,'+-')
   % semilogy(delta_vec2,eigmins,'+-')
    hold on
   % xticks([7*pi/12, 2*pi/3]);
   % xticklabels({'7pi/12', '2\pi/3'});
    xlabel('$\alpha$','interpreter','latex')
    %xticks([pi/2, 2*pi/3]);
    %xticklabels({'\pi/2', '2\pi/3'});
    grid on
    axis tight
    xlabel('$\alpha$','interpreter','latex')
    %xlabel('$\delta_2$','Interpreter','latex')
    ylabel('Smallest eigval','interpreter','latex')
    
    figure(3)
    semilogy(alpha_vec,dist,'+-')
    %semilogy(delta_vec2,dist,'+-')
    hold on
    xlabel('$\alpha$','interpreter','latex')
    %xlabel('$\delta_2$','Interpreter','latex')
    ylabel('Distance particle 1-3','interpreter','latex')
    grid on
    axis tight
   % xticks([7*pi/12, 2*pi/3]);
   % xticklabels({'3pi/4', '2\pi/3'});
    xlabel('$\alpha$','interpreter','latex')
    %xticks([pi/2, 2*pi/3]);
   % xticklabels({'\pi/2', '2\pi/3'});

    figure(4)
    loglog(dist,eigmins,'+-')
    xlabel('Distance particle 1-3','interpreter','latex')
    ylabel('Smallest eigval','interpreter','latex')
    hold on
    axis tight
    grid on

end

