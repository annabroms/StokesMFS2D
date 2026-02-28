function [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc,debug)
%GETPAIRBASISSTOKES Build pair-basis pseudoinverse factors for 2D Stokes
%pair corrections.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap,Cmap_F] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt)
%   [Uf,Yf,Up,Yp,Cmap,Cmap_F] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc)
%   [Uf,Yf,Up,Yp,Cmap,Cmap_F] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc,debug)
%
% Inputs:
%   q            - P-by-1 complex particle centers. rbase_in_c   - Coarse
%   proxy source nodes for one body at the origin. rbase_in_f   - Fine
%   proxy source nodes for one body at the origin. rimage_pairs - Cell
%   array; rimage_pairs{i,j} are extra enhancing source points
%                  used for pair (i,j).
%   refine       - Cell array; refine{i,j} stores near-contact collocation
%                  parameter values for pair (i,j).
%   pairs        - Npair-by-2 list of close pairs [i, j]. opt          -
%   Options struct. Uses fields such as N_f, a_f, precomp,
%                  N_peanut.
%   project      - Optional logical flag. If true, build additional pair
%                  projection blocks used by mobility variants.
%   Lc           - Optional coarse one-body projector for mobility. debug
%   - Optional logical flag. If true, visualizes pair geometry
%                  (sources and collocation nodes) for each processed pair.
%
% Outputs:
%   Uf, Yf  - Cell arrays with fine pair pseudoinverse factors. For each
%             close pair (i,j), beta_pair = Yf{i,j}*(Uf{i,j}*rhs_pair).
%   Up, Yp  - Cell arrays with peanut-compression pseudoinverse factors
%             (empty when opt.N_peanut == 0 or opt.cmap.
%   Cmap    - Cell array with direct coarse-to-coarse pair maps
%             assembled from the factors above (empty when no peanut map
%             and if opt.cmap is false).
%   Cmap_FU  - Reserved output for coarse-to-force/torque mapping,
%             alternatively coarse-to-RBM mapping (empty when no peanut map
%             and if opt.cmap is false).
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

if nargin < 8 || isempty(project)
    project = 0;
end
if nargin < 9
    Lc = [];
end
if nargin < 10 || isempty(debug)
    debug = false;
end

P = opt.P;  
Uf = cell(P);
Yf = cell(P);
N_peanut = opt.N_peanut; 


%If using peanut compression 
if N_peanut
    if opt.cmap
        Up = [];
        Yp = [];
        Cmap = cell(P);
        Cmap_FU = cell(P);
    else
        Up = cell(P);
        Yp = cell(P);
        Cmap = [];
        Cmap_FU = [];
    end
else
    Up = [];
    Yp = [];
    Cmap = [];
    Cmap_FU = [];
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
            rimage_i = rimage_pairs{i,p2};
            rimage_p2 = rimage_pairs{p2,i};
            rin_pair = [rin_1_f; rimage_i; rin_2_f; rimage_p2];

            if debug
                figure(801);
                clf;
                plot(real(rin_1_f),imag(rin_1_f),'r.','MarkerSize',10);
                hold on;
                plot(real(rin_2_f),imag(rin_2_f),'b.','MarkerSize',10);
                plot(real(q(i)+rout_base),imag(q(i)+rout_base),'ro','MarkerSize',4);
                plot(real(q(p2)+rout_base),imag(q(p2)+rout_base),'bo','MarkerSize',4);
                plot(real(fine_1),imag(fine_1),'r+','MarkerSize',6); 
                plot(real(fine_2),imag(fine_2),'b+','MarkerSize',6);
                if ~isempty(rimage_i)
                    plot(real(rimage_i),imag(rimage_i),'ks','MarkerSize',5);
                end
                if ~isempty(rimage_p2)
                    plot(real(rimage_p2),imag(rimage_p2),'kd','MarkerSize',5);
                end
                plot(real(q(i)),imag(q(i)),'rx','MarkerSize',10,'LineWidth',1.5);
                plot(real(q(p2)),imag(q(p2)),'bx','MarkerSize',10,'LineWidth',1.5);
                axis equal;
                grid on;
                title(sprintf('getPairBasisStokes pair (%d,%d)',i,p2), ...
                    'Interpreter','none');
                drawnow;
            end
            
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
                Lf_pair = [];
            end
            
            %% Compute fine basis pseudoinverse
            %Takes in fine grid with image enhancement and fine collocation
            %points.
         
            [Uf_pair,Yf_pair] = getPairBlockStokes(rin_pair,rout_f,Lf_pair,Lr_pair);

            %build matrix to get the evaluation of coarse Stokeslets on one
            %particle in the pair at the time, zero on the other
            Npair = evaluateCoarseOnPair([q(i),q(p2)],rbase_in_c,rout_f);
            Uf{i,p2} = -Uf_pair'*Npair; 

            Yf{i,p2} = Yf_pair; 

            if N_peanut %If true: compress on peanut described by N_peanut points, i.e. find mapping for equivalent coarse sources
                
                peanut_debug = 0; 
                rout_peanut = createPeanut(q(i),q(p2),N_peanut,peanut_debug);    
                rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];      
                [DC,YC] = getPeanutBlockStokes(rin_pair_c,rin_pair,rout_peanut,Lc_pair,Lf_pair); 
  
                if opt.cmap
                    % Determine coarse to coarse map for the pair
                    Cmap{i,p2} = YC*(DC*Yf_pair*(Uf_pair'*Npair)); 
                    %Construct mapping also for the 
                    % i) fource and torque vector (ONLY for resistance), or
                    % ii) RBM vector (ONLY for mobility)
                    Kft_pair = getKftPair(Kf1,Kf2); 
                    Cmap_FU{i,p2} = Kft_pair*Yf_pair*(Uf_pair'*Npair);                   
                else
                    % Store compression for the fine grid
                    Up{i,p2} = DC;
                    Yp{i,p2} = YC;
                end
 
            end
             
    
        end
        
    
    end
    

end

end


