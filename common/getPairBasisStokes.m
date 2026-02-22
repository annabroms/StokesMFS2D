function [Uf,Yf,Up,Yp,Cmap,Cmap_F] = getPairBasisStokes(q,rads,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc)
%getPairBasisStokes computes pair corrections for a 2-body preconditioned
%Stokes MFS solver for circular particles; it solves two LSQ problems via
%SVDs per identified close particle pair and stores the factorisations for
%a backward stable apply. Fine sources are determined via a fine BVP for each 
%pair, and the equivalent coarse sources are computed via matching on a peanut 
%boundary.
%
%Syntax: [Ubf,Ybf,Ucf,Ycf,Cmap,Cmap_F,nimage] = getPairBasisStokes(q,N_f,a_f,rads,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc,Lf,Kf)
%
%Input: 
% q      - Vector of complex valued center coordinates for the P particles
% N_f    - Number of proxy sources on the fine grid
% a_f    - Upsampling factor for the fine grid so that the number of uniform
%         fine collocation points is M_f = N_f*a_f
% rads   - vector of size P containing the particle radii
% rbase_in_c - Coarse grid of proxy sources on a single particle centered at the origin.
% rbase_in_f - Fine grid of proxy sources on a single particle centered at the origin.
% rimage_pairs - A cell array with cell {i,j} containing the image points for
%             particle i in its near contact with particle j
% refine     - A cell array with cell {i,j} containing the parameter values t 
%             of collocation points on particle i in its near contact with 
%             particle j
% pairs     - Matrix with two columns containing a list of all pairs
%             considered close. This can include more pairs than those in 
%             need of image enhancement. 
% opt       - struct containting parameters
% project   - Boolean, true if solving mobility problem
% Lc        - (optional) Projection matrix for single body coarse sources 
%             onto the space of RBM, only for mobility
% Lf        - (optional) Projection matrix for single body fine sources 
%             onto the space of RBM, only for mobility
% Kf        - (optional) Kf' maps fine sources to force and torque, only for mobility
%
% Output:
%
% Uf        - Cell array containing in cell {i,j} the matrix of left
%             singular vectors for the fine evaluation on pair (i,j) multiplied by the
%             evaluation matrix from the coarse sources on a fine grid
% Yf        - Cell array containing in cell {i,j} the matrix formed by VS^+, with S^⁺ a diagonal matrix containing
%             1/sigma for each singular value sigma that is above a set tol and V the 
%             matrix formed by right singular vectors, computed from the
%             fine evaluation on pair (i,j) (fine sources, fine collocation
%             points. Determines the evaluation of fine sources for the
%             pair togehter with Uf so that for each pair, beta =
%             Yf*(Uf*lambda)
%
% Up        - Cell array containing in cell {i,j} the matrix of left
%             singular vectors for the evaluation of coarse sources on the peanut for pair (i,j) multiplied by the
%             evaluation matrix from the fine sources on the peanut grid
% Yp        - Cell array containing in cell {i,j} the matrix formed by VS^+, with S^⁺ a diagonal matrix containing
%             1/sigma for each singular value sigma that is above a set tol and V the 
%             matrix formed by right singular vectors, computed from the
%             evaluation of coarse sources on the peanut corresponding to
%             pair (i,j). Determines coarse sources equivalent to fine
%             sources, together with Up. For each matrix in Up,
%             lambda_effective = Yp*(Up*beta), with beta the fine sources.
%
% Cmap      - Compresses the four matrix blocks above together so that in
%             cell {ij}, Yp*(Up*Yf*(Uf)), which maps lambda_effective
%             <-lambda_coarse
% Cmap_F    - A similar effective mapping for the contribution to the
%             fource and torque
%
% Note: Peanut compression is done if opt.NPeanut > 0. Only particles with
% unit radius supported so far.
%
% Anna Broms Feb 13, 2026

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

            % image locations
            %rimage = [rimage_pairs{i,p2}; rimage_pairs{p2,i}];
       

            % collocation points for the pair
            nout = ceil(a_f*N_f); 
          % nout = ceil(a2*(opt.N_c+size(rimage,1))); %Old choice
            t = linspace(0,2*pi,nout+1);
            t = t(1:end-1)';
            rout_base = rads(k)*(cos(t)+1i*sin(t));
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
%                 rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,i,p2); 
%                 Nother = singleLayer(rbase_in_c+q(i),rout_fine_other,1);
%                 R2 = -Nother*tau_mapped; %read off on particle 2
%                  
%                 rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,p2,i); 
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
                [DC,YC] = getPeanutBlockStokes(rin_pair_c,rin_pair,rout_peanut,rimage,Lc_pair,Lf_pair); 
                          

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





