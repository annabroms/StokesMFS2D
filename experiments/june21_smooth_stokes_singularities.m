clear;
close all
%We test a suitable combination of fundamental solutions for a slender
%ellipse, discretized with a proxy curve (p) or a line of source points
%(l). The slenderness of the ellipse can be adjusted in
%"june21_smooth_stokes.m". In that function, the b.c. can also be changed.
%Stokeslets or Stokeslets + rotlets don't work well for the line
%discretization. 

smat = dec2bin(0:31) - '0';


plotdomain = 0; 


%N = 100;
%N = 900;
%N = N*3/4;
srcvec = [2 1 2 2 4];


for k = 1:size(smat,1) %Loop over combinations of fundamental solutions
    k

    % Make sure the total number of unknowns is approximately the same
    tot_src = sum(smat(k,:).*srcvec);
    N = ceil(900/tot_src);

    if (k<17) || (smat(k,3) && smat(k,5))
        perr(k) = inf;
        lerr(k) = inf;
        solve_perr(k) = inf;
        solve_lerr(k) = inf;
        coeff_p(k) = inf;
        coeff_l(k) = inf;
        continue;
    end

    [err_p,err_l,reserr_p,reserr_l,coefnorm_p,coefnorm_l] = june21_smooth_stokes(smat(k,:),N,plotdomain);
    perr(k) = err_p;
    lerr(k) = err_l;
    solve_perr(k) = reserr_p;
    solve_lerr(k) = reserr_l;
    coeff_p(k) = coefnorm_p;
    coeff_l(k) = coefnorm_l;

end
%%
figure()
subplot(1,2,1)
semilogy(perr,'r.','MarkerSize',15); %visuals surface residual in new nodes
hold on
semilogy(lerr,'k+','MarkerSize',15);
%semilogy(solve_lerr,'k+','MarkerSize',10); %residual at solve stage
%semilogy(solve_perr,'b+','MarkerSize',10);
grid on

subplot(1,2,2)
semilogy(coeff_p,'r.','MarkerSize',15);
hold on
semilogy(coeff_l,'k+','MarkerSize',15);
grid on