function [Uf,Yf,Up,Yp,Cmap,Cmap_F] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc)
%GETPAIRBASISSTOKES Build pair-basis pseudoinverse factors for 2D Stokes pair corrections.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap,Cmap_F] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt)
%   [Uf,Yf,Up,Yp,Cmap,Cmap_F] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc)
%
% Inputs:
%   q            - P-by-1 complex particle centers.
%   rbase_in_c   - Coarse proxy source nodes for one body at the origin.
%   rbase_in_f   - Fine proxy source nodes for one body at the origin.
%   rimage_pairs - Cell array; rimage_pairs{i,j} are extra/image points
%                  used for pair (i,j).
%   refine       - Cell array; refine{i,j} stores near-contact collocation
%                  parameter values for pair (i,j).
%   pairs        - Npair-by-2 list of close pairs [i, j].
%   opt          - Options struct. Uses fields such as N_f, a_f, precomp,
%                  N_peanut.
%   project      - Optional logical flag. If true, build additional pair
%                  projection blocks used by mobility variants.
%   Lc           - Optional coarse one-body projector for mobility.
%
% Outputs:
%   Uf, Yf  - Cell arrays with fine pair pseudoinverse factors. For each
%             close pair (i,j), beta_pair = Yf{i,j}*(Uf{i,j}*rhs_pair).
%   Up, Yp  - Cell arrays with peanut-compression pseudoinverse factors
%             (empty when opt.N_peanut == 0).
%   Cmap    - Cell array with direct coarse-to-coarse pair maps
%             assembled from the factors above (empty when no peanut map).
%   Cmap_F  - Reserved output for coarse-to-force/torque mapping; currently
%             left empty unless provided by downstream extensions.
%
% Notes:
%   - Pair source ordering follows the assembled pair vectors in this file:
%     particle i fine sources, i->j extra sources, particle j fine
%     sources, j->i extra sources, with x-components stacked before y.
%   - If opt.precomp is true, Uf{i,j} includes the coarse-to-fine
%     evaluation operator (Npair) so runtime application avoids explicitly
%     rebuilding that dense block.
%
% Anna Broms, Feb 13, 2026

if nargin<9
    Lc = [];
    project = 0;
end



P = length(q); 
Uf = cell(P);
Yf = cell(P);
N_peanut = opt.N_peanut; 


%If using peanut compression 
if N_peanut
    Up = cell(P);
    Yp = cell(P);
    Cmap = cell(P);
    Cmap_F = cell(P);
else
    Up = [];
    Yp = [];
    Cmap = [];
    Cmap_F = [];
end


%projection matrices the same for every pair as the same coarse and fine
%grid is used for everybody. There are four blocks in L, xx xy, yx, yy...
%Need to build projection for the pair the same way

if ~isempty(Lc)
    Lc_pair = getILpair(Lc);
else
    Lc_pair = [];
    Lf_pair = [];
end

a_f = opt.a_f;
N_f = opt.N_f; 

for i = 1:P

    if ~isempty(pairs)
        neigh = find(pairs(:,1)==i);
        rin_1_f = q(i)+rbase_in_f; %fine grid on first particle in pair

        for k = 1:length(neigh)
            %% get fine grid for a pair 
            %fine grid on second particle in pair
            p2 = pairs(neigh(k),2);  
            rin_2_f = q(p2)+rbase_in_f;        

            % collocation points for the pair
            nout = ceil(a_f*N_f); 
          % nout = ceil(a2*(opt.N_c+size(rimage,1))); %Old choice
            t = linspace(0,2*pi,nout+1);
            t = t(1:end-1)';
            rout_base = cos(t)+1i*sin(t);
            fine_1 = refine{i,p2};
            fine_2 = refine{p2,i};
            rout_f = [q(i)+rout_base; fine_1 ;q(p2)+rout_base; fine_2];
            
            %fine grid of Stokeslets
            rin_pair = [rin_1_f; rimage_pairs{i,p2}; rin_2_f; rimage_pairs{p2,i}];
            
            %% Projection trick (Mobility only)
            % Need the matrix that maps fine sources to rigid body
            % velocities to close the system 
            if project
                B1 = getKmat2D([q(i)+rout_base; fine_1],q(i));
                B2 = getKmat2D([q(p2)+rout_base; fine_2],q(p2));
                Kf1 = getKmat2D(rin_pair(1:end/2),q(i));
                Kf2 = getKmat2D(rin_pair(end/2+1:end),q(p2));
                Lr_pair = getLrPair(B1,B2,Kf1,Kf2);  
                Lf_pair = getLfPair(Kf1,Kf2); 
            else
                Lr_pair = []; 
            end
            
            %% Compute fine basis pseudoinverse
            %Takes in fine grid with image enhancement and fine collocation
            %points. 
            [Uf_pair,Yf_pair] = getPairBlockStokes(rin_pair,rout_f,Lf_pair,Lr_pair);

            if opt.precomp

                %build matrix to get the evaluation of coarse Stokeslets on one
                %particle in the pair at the time, zero on the other
                Npair = evaluateCoarseOnPair([q(i),q(p2)],rbase_in_c,rout_f);
 
%                 % DEBUG
%                 mapped = rand(opt.N_c*2,1);
%                 tau_mapped= rand(opt.N_c*2,1);
%                 rout_fine_other = getFineOther(opt.a_f,opt.N_f,refine,q,i,p2); 
%                 Nother = singleLayer(rbase_in_c+q(i),rout_fine_other,1);
%                 R2 = -Nother*tau_mapped; %read off on particle 2
%                  
%                 rout_fine_other = getFineOther(opt.a_f,opt.N_f,,refine,q,p2,i); 
%                 Nother2 = singleLayer(rbase_in_c+q(p2),rout_fine_other,1);
%                 R1 = -Nother2*mapped; %read off on particle 1
%               
%                 rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)];
%                 rhs2 = -Npair*[tau_mapped(1:end/2); mapped(1:end/2); tau_mapped(end/2+1:end); mapped(end/2+1:end)];

                %rhs and rhs2 should be the same!

                Uf{i,p2} = -Uf_pair'*Npair; 
            else
                Uf{i,p2} = Uf_pair';
            end

            Yf{i,p2} = Yf_pair; 

            if N_peanut %If true: compress on peanut described by N_peanut points, i.e. find mapping for equivalent coarse sources
                
                debug = 0; 
                rout_peanut = createPeanut(q(i),q(p2),N_peanut,debug);    
                rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];
               % [DC,YC] = getPeanutBlock(rin_pair_c,rin_pair,rout_peanut,[nimage{i,p2}; nimage{p2,i}],rimage,s,Lc_pair,Lf_pair); 
                [DC,YC] = getPeanutBlockStokes(rin_pair_c,rin_pair,rout_peanut,Lc_pair,Lf_pair); 
                          

                Up{i,p2} = DC;
                Yp{i,p2} = YC;
                % Determine coarse to coarse map for the pair
                Cmap{i,p2} = -YC*(DC*Yf_pair*(Uf_pair'*Npair)); 

                %Construct mapping also to the fource and torque vector
                %Cmap_F{i,p2} = 


            end
             
    
        end
        
    
    end
    

end

end




