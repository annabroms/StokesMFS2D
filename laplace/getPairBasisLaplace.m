function [Uf,Yf,Up,Yp,Cmap] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt)
%GETPAIRBASISLAPLACE Build pair-basis pseudoinverse factors for scalar Laplace.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt)
%
% Anna Broms, Mar 2026

P = opt.P;
N_f = opt.N_f;
a_f = opt.a_f;
N_peanut = opt.N_peanut;
if ~isfield(opt,'alpha') || isempty(opt.alpha)
    error('getPairBasisLaplace:MissingAlpha','opt.alpha is required.');
end
R = opt.alpha;

if isfield(opt,'show_counter') && ~isempty(opt.show_counter)
    show_counter = logical(opt.show_counter);
else
    show_counter = false;
end

Uf = cell(P);
Yf = cell(P);

if N_peanut
    Up = cell(P);
    Yp = cell(P);
    Cmap = cell(P);
else
    Up = [];
    Yp = [];
    Cmap = [];
end

processed_pairs = 0;
total_pairs = size(pairs,1);

for ii = 1:total_pairs
    i = pairs(ii,1);
    p2 = pairs(ii,2);

    nout = ceil(a_f*N_f);
    t = linspace(0,2*pi,nout+1)';
    t = t(1:end-1);
    rout_base = R*(cos(t)+1i*sin(t));

    fine_1 = refine{i,p2};
    fine_2 = refine{p2,i};
    rout_f = [q(i)+rout_base; fine_1; q(p2)+rout_base; fine_2];

    rimage_i = rimage_vec{i,p2};
    rimage_p2 = rimage_vec{p2,i};
    rin_pair_f = [q(i)+rbase_in_f; rimage_i; q(p2)+rbase_in_f; rimage_p2];

    [Uf_pair,Yf_pair] = getPairBlockLaplace(rin_pair_f,rout_f);

    Npair = evaluateCoarseOnPairLaplace([q(i);q(p2)],rbase_in_c,rout_f);
    Uf{i,p2} = -Uf_pair'*Npair;
    Yf{i,p2} = Yf_pair;

    if N_peanut
        rout_peanut = createPeanut(q(i),q(p2),N_peanut,0,R);
        rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];
        [DC,YC] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut);

        if isfield(opt,'cmap') && opt.cmap
            Cmap{i,p2} = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
        else
            Up{i,p2} = DC;
            Yp{i,p2} = YC;
        end
    end

    processed_pairs = processed_pairs + 1;
    if show_counter
        fprintf('getPairBasisLaplace: processed pair %d/%d (%d,%d)\n', ...
            processed_pairs,total_pairs,i,p2);
    end
end

end
