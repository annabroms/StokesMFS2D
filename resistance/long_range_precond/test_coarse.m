%test coarse resistance
clear;
close all; 
P = 100;
delta = 1; 
N = 60;
M = N*1.2; 
for k = 1:2*P
    blocksN{k} = ones(N,1);
    blocksM{k} = ones(M,1);
end
AN = blkdiag(blocksN{:});
AM = blkdiag(blocksM{:});

opt = get2Dparams(); 
[q,B] = grow_cluster(P,delta,2);
rads = ones(P,1); 
[rout,~,rin,rimage,nimage,pair_points] = get2DImageGrid(q,rads,opt);

G = stokSLPmat(rin,rout,1);

Rc = AM'*G*AN;

