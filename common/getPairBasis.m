function [Uf,Yf,Up,Yp,Cmap,Cmap_F,nimage] = getPairBasis(q,N_f,a_f,rads,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc,Lf,Kf,debug)
%getPairCorrection computes pair corrections for a 2-body preconditioned
%Stokes MFS solver for circular particles; it solves two LSQ problems via
%SVDs per identified close particle pair and stores the factorisations for
%a backward stable apply. Fine sources are determined via a fine BVP for each 
%pair, and the equivalent coarse sources are computed via matching on a peanut 
%boundary.
%
%Syntax: [Ubf,Ybf,Ucf,Ycf,Cmap,Cmap_F,nimage] = getPairBasis(q,N_f,a_f,rads,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,project,Lc,Lf,Kf,debug)
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
% debug     - (optional) Logical flag. If true, visualizes pair geometry
%             (sources and collocation nodes) for each processed pair.
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
% Anna Broms April 4, 2025


if nargin < 11 || isempty(project)
    project = 0;
end
if nargin < 12
    Lc = [];
end
if nargin < 13
    Lf = [];
end
if nargin < 14
    Kf = [];
end
if nargin < 15 || isempty(debug)
    debug = false;
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

s = opt.s;


nimage = cell(P);

%projection matrices the same for every pair as the same coarse and fine
%grid is used for everybody. There are four blocks in L, xx xy, yx, yy...
%Need to build projection for the pair the same way

if ~isempty(Lc)
    Lc_pair = getILpair(Lc);
    Lf_pair = getILpair(Lf);
else
    Lc_pair = [];
    Lf_pair = [];
end


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
            rimage = [rimage_pairs{i,p2}; rimage_pairs{p2,i}];
            ntemp =  randn(size(rimage,1),2); %generate here? Or store from get2Dgrid?
            nimage{i,p2} = (ntemp(1:end/2,1)+1i*ntemp(1:end/2,2))./abs(ntemp(1:end/2,1)+1i*ntemp(1:end/2,2));
            nimage{p2,i} = (ntemp(end/2+1:end,1)+1i*ntemp(end/2+1:end,2))./abs(ntemp(end/2+1:end,1)+1i*ntemp(end/2+1:end,2));

            % collocation points for the pair
            nout = ceil(a_f*N_f); 
          % nout = ceil(a2*(opt.N_c+size(rimage,1))); %Old choice
            t = linspace(0,2*pi,nout+1);
            t = t(1:end-1)';
            rout_base = rads(k)*(cos(t)+1i*sin(t));
            t1 = refine{i,p2};
            fine_1 = q(i)+rads(k)*(cos(t1)+1i*sin(t1));
            t2 = refine{p2,i};
            fine_2 = q(p2)+rads(k)*(cos(t2)+1i*sin(t2));
            rout_f = [q(i)+rout_base; fine_1 ;q(p2)+rout_base; fine_2];
            
            %fine grid of Stokeslets
            rin_pair = [rin_1_f; rin_2_f];

            if debug
                figure(800);
                clf;
                plot(real(rin_1_f),imag(rin_1_f),'r.','MarkerSize',10);
                hold on;
                plot(real(rin_2_f),imag(rin_2_f),'b.','MarkerSize',10);
                plot(real(q(i)+rout_base),imag(q(i)+rout_base),'ro','MarkerSize',4);
                plot(real(q(p2)+rout_base),imag(q(p2)+rout_base),'bo','MarkerSize',4);
                plot(real(fine_1),imag(fine_1),'r+','MarkerSize',6);
                plot(real(fine_2),imag(fine_2),'b+','MarkerSize',6);
                if ~isempty(rimage)
                    plot(real(rimage),imag(rimage),'ks','MarkerSize',5);
                end
                plot(real(q(i)),imag(q(i)),'rx','MarkerSize',10,'LineWidth',1.5);
                plot(real(q(p2)),imag(q(p2)),'bx','MarkerSize',10,'LineWidth',1.5);
                axis equal;
                grid on;
                title(sprintf('getPairBasis pair (%d,%d)',i,p2), ...
                    'Interpreter','none');
            end
            
            %% Projection trick (Mobility only)
            % Need the matrix that maps fine sources to rigid body
            % velocities to close the system 
            if project
                B1 = getKmat2D([q(i)+rout_base; fine_1],q(i));
                B2 = getKmat2D([q(p2)+rout_base; fine_2],q(p2));
                Lr_pair = getLrPair(B1,B2,Kf,Kf);            
            else
                Lr_pair = []; 
            end
            
            %% Compute fine basis pseudoinverse
            %Takes in fine grid with image enhancement and fine collocation
            %points. 
            [Uf_pair,Yf_pair] = getPairBlock([q(i);q(p2)],rin_pair,rout_f,rimage,[nimage{i,p2}; nimage{p2,i}],s,Lf_pair,Lr_pair,opt.proj_all);

            if opt.precomp

                %build matrix to get the evaluation of coarse Stokeslets on one
                %particle in the pair at the time, zero on the other
                Npair = evaluateCoarseOnPair([q(i),q(p2)],rbase_in_c,rout_f);
 
                Uf{i,p2} = -Uf_pair'*Npair; 
            else
                Uf{i,p2} = Uf_pair';
            end

            Yf{i,p2} = Yf_pair; 

            if N_peanut %If true: compress on peanut described by N_peanut points, i.e. find mapping for equivalent coarse sources
                
                peanut_debug = 0; 
                rout_peanut = createPeanut(q(i),q(p2),N_peanut,peanut_debug);    
                rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];
                [DC,YC] = getPeanutBlock(rin_pair_c,rin_pair,rout_peanut,[nimage{i,p2}; nimage{p2,i}],rimage,s,Lc_pair,Lf_pair); 
                          

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




