clear
close all

smat = dec2bin(0:31) - '0';
plotdomain = 0; 
clust_num = 1; 

for i = 1:size(smat,1)

    i
    
    s = smat(i,:); 
    if s(:,3) && s(:,5)
        err(i) = inf;
        reserr(i) = inf; 
        lambda_norm(i) = inf; 

        continue;
    end

    [err_ext,reserr_ext,coefnorm_ext] = june18_polygon_stokes(s,plotdomain);
    err(i) = err_ext;
    reserr(i) = reserr_ext;
    lambda_norm(i) = coefnorm_ext;
    
end

%% Visualise

col = lines(32); 
msize= 15; 

figure()
subplot(1,2,1)
for k = 1:32
    semilogy((k-1)*ones(1,clust_num),err(k),'.','Color',col(k,:),'MarkerSize',msize);
    hold on
    semilogy((k-1)*ones(1,clust_num),reserr(k),'+','Color',col(k,:),'MarkerSize',msize);
end
xlim([0,30])
ylabel('Residual in collocation points \& new points','interpreter','latex')
xlabel('Source types at image points')
sgtitle('Non-convex polygon -- 150 clustered nodes','interpreter','latex')
max(err)
grid on

subplot(1,2,2)
for k = 1:32
    semilogy((k-1)*ones(1,clust_num),lambda_norm(k),'.','Color',col(k,:),'MarkerSize',msize);
    hold on
end
xlim([0,30])
ylabel('$\|\vec \lambda\|$','interpreter','latex')
xlabel('Source types at image points','interpreter','latex')
grid on
max(lambda_norm)