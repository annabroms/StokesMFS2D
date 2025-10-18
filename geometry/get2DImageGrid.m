function [rout, weights, rin, rimage, nimage, pair_points, pairs, rimage_pairs, refine, rin_base] = get2DImageGrid(q, rads, opt)
%GET2DIMAGEGRID Distributes source, collocation, and image points for circular particles in 2D
%
% Syntax:
%   [rout, rin, rimage, nimage, pair_points, pairs, rimage_pairs, refine, rin_base] = get2DImageGrid(q, rads, opt)
%
% Inputs:
%   q          - Complex vector of particle center coordinates (length P)
%   rads       - Vector of particle radii (length P)
%   opt        - Struct with options including:
%                   * Rp_c  - Proxy radius for the coarse grid (the only grid used for 1-body preconditioning)
%                   * N_c   - Number of proxy points on the coarse grid
%                   * a_c   - Upsampling factor for the coarse grid so that the number of coarse collocation points 
%                             is M_c = a_c × N_c
%                   * Rp_f - Proxy radius for the fine grid (used for 2-body preconditioning)
%                   * N_f  - Number of proxy points on the fine grid
%                   * a_f  - Upsampling factor for fine grid (M_f = a_f × N_f)
%                   * pc   - Logical flag for pair corrections: if true,
%                            return only the basic coarse collocation grid as rout
%                   * delta_pair - Proximity threshold for determining particle pairs requiring close interaction pair corrections
%                   * image      - Logical flag: activate image systems for close interactions
%                   * s    - source types at image points, needed here to
%                            determine the number of extra collocation points for
%                            image lines
%                   * M_image - scalar that will set the number of extra
%                               collocation points per close neighbour     
% Outputs:
%   rout         - Complex vector of collocation points on all particles:
%                  body1, body2, ... , body P
%   rin          - Complex vector of source points on all particles, body1, body2, ... , body P
%   rimage       - Complex vector of image source locations for all
%                  particles
%   nimage       - Complex vector of stresslet directions at image points, stored as x + iy
%   pair_points  - P×2 matrix; row j gives total number of image and collocation points for particle j
%   pairs        - Matrix listing particle pairs considered close; may include more than image-enhanced pairs
%   rimage_pairs - Cell array {i,j} of image points for particle i near particle j
%   refine       - Cell array {i,j} of parameter values (on [0,2π)) for collocation points on particle i near j
%   rin_base     - Vector of proxy source points for a reference particle centered at the origin (fine grid)
%
% Description:
%   Constructs MFS discretization geometry for a system of circular particles.
%   Handles both 1-body and 2-body preconditioned setups, and includes image system enhancement
%   for improved accuracy for near close-contact interactions. 
%
% Notes:
%   - For 1-body preconditioning, we set a_f = a_c and N_f = N_c
%   - For 2-body preconditioning (pair corrections), set pc = true.
%     In this case, image-related outputs (rimage, nimage) will not be used.
%   - The "image points" are approximations and discretise the line of images from the
%     image accumulation point to the fine proxy radius for each interaction with
%     half a Chebyshev grid clustered towards the accumulation point.
%   - A graded periodic trapezoid rule is used for extra collocation points
%     with clustered nodes above any image line.
%
% Anna Broms, April 9, 2025
%

% Unpack parameters
Rp_c = opt.Rp_c;
N_c = opt.N_c;
a_c = opt.a_c;
Rp_f = opt.Rp_f;
N_f = opt.N_f;
a_f = opt.a_f; 
image = opt.image;
pc = opt.pc; 
delta_pair = opt.delta_pair;
M_image = opt.M_image; 
if opt.image
    im_types = length(opt.s); %the number of source types at the image points will also determine the number of extra collocation points
end

if ~pc
    N_f = N_c;
    a_f = a_c;
    Rp_f = Rp_c;
end


P = length(q); %Number of particles

pairs = [];
rimage_pairs = [];
refine = [];


if image
    %parameters to give number of image points

    %accumulation point, given Rg and delta. Closed formula from fixed point of reflection formula
    accstop_fine = (1-Rp_f)^2/Rp_f; 
    %accstop_coarse = (1-Rp_c)^2/Rp_c; 
    
    %Empirical linear relationship in log(delta) for the number of images
    %per pair
    slope = (-30+1)/(3+log10(accstop_fine)); %35
    
    %slope = (-30+1)/(3+log10(accstop));
    m = 1-log10(accstop_fine)*slope;
    
    %Set local upsampling on the boundary (empirical)
    a_loc = 6; 
    %a_loc = 12; 
    a_loc2 = 60;
    a_loc2 = 30; 
    %a_loc2 = 6;

    %a_loc = 1;
    
    %store image points per pairs, refined collocation points per pair
    rimage_pairs = cell(P,P);
    refine = cell(P,P);
    
    %accstop_fine = 1; %debug; 
  
    for i = 1:P
        ind = setdiff(1:P,i); %check neighbours with everyone else
        for k = ind
            delta = abs(q(i)-q(k))-rads(k)-rads(i);
            d = delta/rads(i);
            %if delta< Rp*rads(i)
           
            if delta < accstop_fine*rads(i)
                %new close pair detected
                if i<k
                    pairs = [pairs; i k]; 
                end
                %% Assign image points    
                if isfield(opt, 'n_clusters')
                    sample_nbr = opt.n_clusters; 
                else
                    %get number of image locations to be generated, from empiric
                    %relationship
                    sample_nbr = max(ceil(slope*log10(d)+m),0); 
                end

                if ~sample_nbr
                    continue;
                end
        
                %generate points from Rp out to the accumulation point.
                %first determine accumulation point            

                a1  = rads(i); 
                a2 = rads(k); 
                D = a1+a2+delta;
                S1 = D^2-a2^2+a1^2; 
                xa1 = (S1/D/2-sqrt((S1^2)/D^2/4-a1^2))/D;
                acc = xa1*D; %accumlation point

                %Compute half a Chebyshev grid clustered towards the
                %accumlation point
                n = 2*sample_nbr;  
                xx = cos(pi*(2*(0:n)'+1)./(2*n+2)); 
                xx = -xx(sample_nbr+1:end);

                %Different discretisation of image line
%                 beta = log(5/sqrt(delta));
%                 %beta = log(3/delta);
%                 [xx2, ~] = gradedptr(n, beta);
%                 xx2 = real(xx2)/(pi);
%                 xx = xx2(end/2+2:end)'-1;

                
                
                aa = rads(i)*Rp_f*1.01; %start just exterior to the proxy grid     
                aa = rads(i)*Rp_f*1.05;
                %aa = 0.6; %test
                
                %Rescale image line to be exactly between the end points
                % (aa, accumulation point)

               % acc = acc+0.01; % Only for debugging purposes!
               % acc = 0.99;

                %different sampling
                xx = logspace(2*log10(sqrt(delta)),0,sample_nbr)';

                %different sampling
                % s = sqrt(1:sample_nbr) - sqrt(sample_nbr);
                % A_clust = 0.6;
                % sigma = 4;
                % xx = A_clust * exp(sigma*s)';

                xx = xx.*(acc-aa);
                %t = ((xx - 0) * (acc-aa) / (max(xx) - 0)) + aa;
                t = acc-xx+xx(1);


                %t = ((xx - 0) * (acc - aa) / (max(xx) - 0)) + aa;


                %xx = aa*logspace(2*log10(sqrt(delta)),0,sample_nbr)';
                %t = [acc; acc-xx]; %large resiudal

                mid = (acc+aa)/2;
                
             
                % If assigning two lines of image points, meeting at an
                % angle alphas
                alpha = opt.alpha; %pi-alpha;
                t_new = t*cos(alpha)+1i*sin(alpha)*t+mid;
                t_new = t_new-t_new(1)+t(1);
                %figure()
                %plot(real(t_new),imag(t_new),'+')


                %t = acc; %debug

                % Assign image poins on particle i, close to touching
                % particle k
                line = q(i)+t*(q(k)-q(i))./abs(q(k)-q(i));
                %line1 = q(i)+t_new*(q(k)-q(i))./abs(q(k)-q(i));
                %line2 = q(i)+conj(t_new)*(q(k)-q(i))./abs(q(k)-q(i));


               % line2 = []; 
               % hold on
               % plot(real(line1),imag(line1),'*')
               % plot(real(line2),imag(line2),'*')

               %% Test to put the extra sources along an arc instead.
               %t = logspace(-6,log10(alpha),sample_nbr)';
              % t = [-t; 0; t];
               %line = q(i)+acc*(cos(t)+1i*sin(t))*(q(k)-q(i))./abs(q(k)-q(i));

               %% Store
                
                 % figure()
                 % plot(real(line),imag(line),'*');
                 % axis equal
                rimage_pairs{i,k} = line;
                %rimage_pairs{i,k} = [line1; line2];



        
                %% Assign collocation points
                %Need to store also the number of extra collocation  points
    
                %location for the point of contact. 
                xstar = q(i)+rads(i)*(q(k)-q(i))./abs(q(k)-q(i));

                %determine the corresponding angle (parameter
                %value) at point of contact
                tstar = atan2(imag(xstar-q(i)), real(xstar-q(i))); 

                
                %Tried a bunch of different strategies for assigning the
                %extra collocation points... 
                % * Two equally spaced intervals with different densities
                % * Logaritmically spaced points towards tstar...

                  %t= (tstar-h*sample_nbr):h:(tstar+h*sample_nbr);
        
                  %Set points with a certain density... 
%                 nout = N_f+2*2*(sample_nbr+1);
%                 %nout = 40+2*2*(sample_nbr+1); 
%                 h = 2*pi/nout/a_loc; 
%                 h2 = 2*pi/nout/a_loc2; 
%                 
%                 %the extra sample points on particle k
%                 t2 = (tstar-5*h2*sample_nbr):h2:(tstar+5*h2*sample_nbr);
%                 t= (tstar-4*h*sample_nbr):h:(tstar+4*h*sample_nbr); %seems to work well for Nc = 40... 
% 
%                 %More clustered... 
%                 strip = logspace(-6,log10(4*h*sample_nbr),10*sample_nbr);
%                 %strip = logspace(-6,log10(25*h*sample_nbr),sample_nbr);
%                 %strip = logspace(-6,log10(25*h*sample_nbr),10*sample_nbr);
%                 t = [tstar-strip tstar+strip];
%                 t = sort(t);

                %% Instead, use a graded PTR (see Alex's BIE book, Chapter 4)
               % beta = log(3/sqrt(delta));
                %beta = log(3/delta);
                beta = log(5/sqrt(delta));
                N = M_image*round(beta);
                N = max([N,ceil(1.2*im_types*sample_nbr)]);
                [t, ~] = gradedptr(N, beta);

                %test alternative
                M_clust = sample_nbr*4*sum(opt.s);
                t = 0.5*[logspace(-6,0,M_clust) -logspace(-6,0,M_clust)];
                %must be rotated
                t = t+tstar;

                %Store the additional collocation points for this pair
                refine{i,k} = t';                
               
            elseif pc %if pair-corrections  
                if (delta < delta_pair)
                    %want to use fine grid of stokeslets for these distances, 
                    % but no image points as the particles are not
                    % sufficiently close
                    %and only a uniform grid at the boundary.
                    
                    if i<k
                        pairs = [pairs; i k];
      
                    end
                end
            end
    
        end
    end
end

%% create inner proxy grid
tin = linspace(0,2*pi,N_f+1);
tin = tin(1:end-1)';
tin = tin+0.2;
rin = [];
rimage = [];
nimage = [];
rout = []; 
weights = [];

Mf = ceil(a_f*N_f); 

for k = 1:P
    %Assign inner grid to particles, based on fine grid (if 1-body precond,
    %same as coarse grid)
    rin_base = Rp_f*rads(k)*cos(tin)+1i*Rp_f*rads(k)*sin(tin);
    rin = [rin; q(k)+rin_base];

    if image
        %start with the fine uniform grid
        t = linspace(0,2*pi,Mf+1);
        t = t(1:end-1)';

        %store global list of all image locations
        % -- only needed with one body precond! 
        rtemp = [];
        for i = 1:P
            rki = rimage_pairs{k,i};
            rtemp = [rtemp; rki];           
        end
        
        %also, construct "normals", used for stresslets 
        num_im = length(rtemp);
        pair_points(k,1) = num_im;

        if ~pc
            ntemp = randn(num_im,1)+1i*randn(num_im,1);
            rimage = [rimage; rtemp];
            nimage = [nimage; ntemp./abs(ntemp)];
        end
        
        %assign outer grids to particles
        %constructed from a base grid plus local refinements. 

        % Add the extra target nodes for particles with near contacts
        t_extra = [];
        for i = 1:P
            tki = refine{k,i}; 
            t_extra = [t_extra; tki];           
        end
        
        if ~pc %if no pair corrections (fine grid to be returned)
            t = [t; t_extra];    
            t = mod(t,2*pi);
            t = sort(t);
            w_k = [diff(t); abs(t(end)-(t(1)+2*pi))]; %weights to be used in left preconditioning
        else 
            %return only the uniform discretization
            t = linspace(0,2*pi,a_c*N_c+1);
            t = t(1:end-1)';
        end
        
        pair_points(k,2) = Mf+length(t_extra); 
 
    else        
        t = linspace(0,2*pi,Mf+1);
        t = t(1:end-1)';

        % uniform grid of collocation points 
        pair_points(k,1) = 0; %no images
        pair_points(k,2) = Mf; %standard number of collocation points
    end
    
    
    rout_part = q(k)+rads(k)*(cos(t)+1i*sin(t));
    rout = [rout; rout_part]; %add to global list
   % weights = [weights; w_k];

end

end