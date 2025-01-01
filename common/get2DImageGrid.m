function [rout,rin,rimage,nimage,pair_points,pairs,rimage_vec,refine,rin_base] = get2DImageGrid(q,rads,Rg_c,a_c,N_c,image,Rg_f,a_f,N_f,basic,delta_pair)
%GET2DIMAGEGRID(q,rads,Rg_c,a_c, N_c, image,Rg_f,a_f,N_f,basic) returns collocation and source points,
%along with a list of image points for each particle. nimage contains
%normals at the image points and pair_points stores data about the number of image
% points and total collocation points for each body. Circles of radii rads
% are centered at q. Rg_c, Rg_f specifies the radius of the proxy circle for 
% the coarse and fine grid (if used), np is the
% number of proxy sources, a the upsampling factor for the uniform grid of collocation 
% points, and "image" specifies if image points are to be used or not.
% basic is a flag to return only the basic collocation points on all bodies
% if true, otherwise the refined grid. 


if nargin <7
    Rg_f = Rg_c; 
    a_f = a_c; 
    N_f = N_c; 
    basic = 0; 
    pc = 0;
elseif nargin < 8
    pc = 0; %no pair-corrections
else
    pc = 1; 
end

P = length(q); 

if image
    %parameters to give number of image points

    %accumulation point, given Rg and delta. Closed formula from fixed point of reflection formula
    accstop_fine = (1-Rg_f)^2/Rg_f; 
    %accstop_coarse = (1-Rg_c)^2/Rg_c; 
    
    %Empirical linear relationship in log(delta) for the number of images
    %per pair
    slope = (-35+1)/(3+log10(accstop_fine));
    
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
    rimage_vec = cell(P,P);
    refine = cell(P,P);
    pairs = [];
    %Assign image points
    for i = 1:P
        ind = setdiff(1:length(q),i);
        for k = ind
            delta = abs(q(i)-q(k))-rads(k)-rads(i);
            d = delta/rads(i);
            %if delta< Rg*rads(i)
            if delta < accstop_fine*rads(i)
                if i<k
                    pairs = [pairs; i k];
                end
                %get number of points to be generated
                sample_nbr = ceil(slope*log10(d)+m); 
                %sample_nbr = 30; 
                %sample_nbr = ceil(slope*log10(d)+m);  %just for debugging
    
                %test:
                %sample_nbr = sample_nbr+10; 
        
                %generate points from Rg to the accumulation point.
                %first determine accumulation point
                
                %NB: this formula holds only for equal cirlces
                xa11 = ((1+d/2)-sqrt(d+(d^2)/4))/(2+d);
                
                a1  = rads(i);
                a2 = rads(k); 
                D = a1+a2+delta;
                S1 = D^2-a2^2+a1^2;
                S2 = -D^2+a1^2-a2^2;
                %xa1 = (S1/D/2-sqrt((S1^2)/D^2/4-a1^2))/(a1+a2+d);
                %xa2 = (S2/D/2+sqrt((S2^2)/D^2/4-a2^2))/(a1+a2+d);
                xa1 = (S1/D/2-sqrt((S1^2)/D^2/4-a1^2))/(a1+a2+delta);
                xa2 = (S2/D/2+sqrt((S2^2)/D^2/4-a2^2))/(a1+a2+delta);

                %Compute chebyshev
                n = 2*sample_nbr; 
        
                xx = cos(pi*(2*(0:n)'+1)./(2*n+2)); 
                xx = -xx(sample_nbr+1:end);
                aa = rads(i)*Rg_f*1.01;
                %bb = xa1*(a1+a2+d);
                bb = xa1*(a1+a2+delta);
                %t = 0.5*(a+b)+0.5*(b-a)*xx;                
                t = ((xx - 0) * (bb - aa) / (max(xx) - 0)) + aa;
                %t = t*rads(i);
        
                rimage_vec{i,k} = q(i)+t*(q(k)-q(i))./abs(q(k)-q(i));
        
                %Need to store also the number of extra collocation  points
    
                %location for the point of contact. 
                xstar = q(i)+rads(i)*(q(k)-q(i))./abs(q(k)-q(i));
                %determine the corresponding angle (parameter
                %value) at point of contact
                tstar = atan2(imag(xstar-q(i)), real(xstar-q(i))); 

                %Tried a bunch of different strategies for assigning the
                %extra collocation points... 

                %t= (tstar-h*sample_nbr):h:(tstar+h*sample_nbr);
        
                nout = N_f+2*2*(sample_nbr+1);
                %nout = 40+2*2*(sample_nbr+1); 
                h = 2*pi/nout/a_loc; 
                h2 = 2*pi/nout/a_loc2; 
                
                %the extra sample points on particle k
                %t= (tstar-h*sample_nbr):h:(tstar+h*sample_nbr); %obs doubled number
                t2 = (tstar-5*h2*sample_nbr):h2:(tstar+5*h2*sample_nbr);
               % t2 = [];
                t= (tstar-5*h*sample_nbr):h:(tstar+5*h*sample_nbr); %obs doubled number
                t= (tstar-3*h*sample_nbr):h:(tstar+3*h*sample_nbr);
                t= (tstar-3.5*h*sample_nbr):h:(tstar+3.5*h*sample_nbr); %depends on resolution probably... 
                t= (tstar-4*h*sample_nbr):h:(tstar+4*h*sample_nbr); %seems to work well for Nc = 40... 

                strip = logspace(-6,log10(4*h*sample_nbr),sample_nbr);
                t = [tstar-strip tstar+strip];
                t = sort(t);

                % t= (tstar-25*h*sample_nbr):h:(tstar+25*h*sample_nbr);

                 %strip = logspace(-6,log10(4*h*sample_nbr),sample_nbr);
                 strip = logspace(-6,log10(25*h*sample_nbr),sample_nbr);
                 %t = [tstar-strip tstar+strip];
                 %t = sort(t);
% 
%                 %The following seems to work pretty well:
%                  strip = logspace(-6,log10(4*h*sample_nbr),8*sample_nbr);
%                  strip = logspace(-6,log10(8*h*sample_nbr),8*sample_nbr);
%                  strip = logspace(-6,log10(25*h*sample_nbr),25*sample_nbr);
%                  strip = logspace(-6,log10(4*h*sample_nbr),sample_nbr);
%                  t = [tstar-strip tstar+strip];
%                  t = sort(t);
               % t= (tstar-5*h*sample_nbr):h:(tstar+5*h*sample_nbr); %obs doubled number
%                 if max(t)<0
%                     t = t+2*pi;
%                     t2 = t2+2*pi; 
%                 end
               

                %Store the additional collocation points for this pair
                refine{i,k} = [t'; t2'];
                
               
            elseif pc %if pair-corrections  
                if (delta < delta_pair)
                    % what is counted as a pair (a larger gap)
                    %want to use fine grid of stokeslets for these distances, 
                    % but no image points,
                    %and only a uniform grid at the boundary.
                    
                    if i<k
                        pairs = [pairs; i k];
      
                    end
                end
            end
    
        end
    end
end

%% create inner grid
tin = linspace(0,2*pi,N_f+1);
tin = tin(1:end-1)';
% rin = [Rg*cos(t); Rg*sin(t)]';
rin = [];
rimage = [];
nimage = [];
rout = []; 

for k = 1:P
    %Assign inner grid to particles, based on fine grid
    rin_base = Rg_f*rads(k)*cos(tin)+1i*Rg_f*rads(k)*sin(tin);
    rin = [rin; q(k)+rin_base];
    if image
    %... and images
        rtemp = [];
        for i = 1:P
            rki = rimage_vec{k,i};
            rtemp = [rtemp; rki];           
        end
        rimage = [rimage; rtemp];
    
        %Construct "normals"
        num_im = length(rtemp);
        pair_points(k,1) = num_im;
        ntemp = randn(num_im,1)+1i*randn(num_im,1);
        nimage = [nimage; ntemp./abs(ntemp)];
        
        %assign outer grids to particles
        %constructed from a base grid plus local refinements. 
    
        %nout = ceil(a_f*(N_f+2*num_im)); 
        nout = ceil(a_f*N_f); 
        %nout = ceil(a*(np+num_im)); 
    else
        nout = ceil(a_c*N_f); %fine grid ?? 
        pair_points(k,1) = 0;
        pair_points(k,2) = nout;
        n_part(k) = nout;
    end
    
    t = linspace(0,2*pi,nout+1);
    t = t(1:end-1)';
    %rout_base = rads(k)*(cos(t)+1i*sin(t));
    %rout = [rout; q(k)+rout_base];
    
    if image
    % Add the extra target nodes
        %t = refine{k};
        t_extra = [];
        for i = 1:P
            tki = refine{k,i}; %Some mapping to be done... 
    %         if (max(t_extra))<0
    %             t_extra = t_extra+2*pi; 
    %             refine{k} = t_extra;
    %         end
            t_extra = [t_extra; tki];           
        end
        
        t = [t; t_extra];    
        t = sort(t);
        %tt{k} = t; %for debugging; 

        n_part(k) = nout+ length(t_extra); 
        pair_points(k,2) = n_part(k);

        %rout = [rout; q(k)+rads(k)*(cos(t)+1i*sin(t))];
       
    end
    if basic
        t = linspace(0,2*pi,ceil(a_c*N_c)+1)';
        t = t(1:end-1); 
        rout_part = q(k)+rads(k)*(cos(t)+1i*sin(t));
    else
        rout_part = q(k)+rads(k)*(cos(t)+1i*sin(t));
    end
    rout = [rout; rout_part]; 

end

end