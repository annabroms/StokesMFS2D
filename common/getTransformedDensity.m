function [tau_stokes_x,tau_stokes_y,rot,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(tau,rimage,Uii,Yii,P,N_small,PM,pairs,s,ver)

if nargin<10
    ver = 1; 
end

if isempty(rimage)
    images = 0;
    UU = Uii{1}; %if same pseudo-inverse for all particles
    Y = Yii{1};
else
    images = 1; 
end
 
%assuming stresslets and potential dipoles at image
rot = zeros(sum(pairs(:,1)),1);
tau_stress_x = zeros(sum(pairs(:,1)),1); 
tau_stress_y = zeros(sum(pairs(:,1)),1); 
tau_pot_x = zeros(sum(pairs(:,1)),1); 
tau_pot_y = zeros(sum(pairs(:,1)),1); 

% tau_rot_x = zeros(sum(pairs(:,1)),1); 
% tau_rot_y = zeros(sum(pairs(:,1)),1); 


%map densities back
%First, map density for all particles
if s(1) %if also stokeslets at image points
    tau_stokes_x = zeros(P*N_small+sum(pairs(:,1)),1);
    tau_stokes_y = zeros(P*N_small+sum(pairs(:,1)),1);
else
    tau_stokes_x = zeros(P*N_small,1);
    tau_stokes_y = zeros(P*N_small,1);
end

start_ind = 0;
start_colloc = 0;

if ~ver %determines ordering of unknowns in pseudoinverse
    for i = 1:P
        %Retrieve self evaluation blocks
        if images
            UU = Uii{i};
            Y = Yii{i};
        end
        tau_particle_x = tau(start_colloc+1:start_colloc+pairs(i,2));
        tau_particle_y = tau(PM+start_colloc+1:PM+start_colloc+pairs(i,2));
        step1 = UU*[tau_particle_x;tau_particle_y]; %here x and y follow each other
        tau_mapped = Y*step1; %this is the mapped density for this particle to throw in to the kernel(s)
        
            
        tau_stokes_x((i-1)*N_small+1:i*N_small) = tau_mapped(1:N_small);      
        tau_stokes_y((i-1)*N_small+1:i*N_small) = tau_mapped(N_small+2*pairs(i,1)+1:2*N_small+2*pairs(i,1));

        
        if images
            tau_stress_x(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(N_small+1:N_small+pairs(i,1));
            tau_pot_x(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(N_small+1+pairs(i,1):N_small+2*pairs(i,1));
            
            tau_stress_y(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(2*N_small+2*pairs(i,1)+1:2*N_small+3*pairs(i,1));
            tau_pot_y(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(2*N_small+3*pairs(i,1)+1:end);           
        
            %the number of images vary!
            start_ind = start_ind+pairs(i,1);
            
        end
        start_colloc = start_colloc+pairs(i,2);
    
    end


else
    for i = 1:P
        %Retrieve self evaluation blocks
        if images
            UU = Uii{i};
            Y = Yii{i};
        end
        tau_particle_x = tau(start_colloc+1:start_colloc+pairs(i,2));
        tau_particle_y = tau(PM+start_colloc+1:PM+start_colloc+pairs(i,2));
        step1 = UU*[tau_particle_x;tau_particle_y]; %here x and y follow each other
        tau_mapped = Y*step1; %this is the mapped density for this particle to throw in to the kernel(s)
    
        %Split, depending on which fundamental solution
        if ~s(1)
            
            tau_stokes_x((i-1)*N_small+1:i*N_small) = tau_mapped(1:N_small);
            tau_stokes_y((i-1)*N_small+1:i*N_small) = tau_mapped(N_small+1:2*N_small);
            
            if images
                tau_stress_x(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(2*N_small+1:2*N_small+pairs(i,1));
                tau_stress_y(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(2*N_small+pairs(i,1)+1:2*N_small+2*pairs(i,1));
                tau_pot_x(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(2*N_small+1+2*pairs(i,1):2*N_small+3*pairs(i,1));
                tau_pot_y(start_ind+1:start_ind+pairs(i,1)) = tau_mapped(2*N_small+3*pairs(i,1)+1:end);           
            
                %the number of images vary!
                start_ind = start_ind+pairs(i,1);
                
            end
        else
            
            warning('To be written!')
            tau_stokes_x((start_ind+1:start_ind+pairs(i,1))) = tau_mapped(1:N_small);
            tau_stokes_y((i-1)*N_small+1:i*N_small) = tau_mapped(N_small+1:2*N_small);
    
        
        
    
    
        
        
        end
        start_colloc = start_colloc+pairs(i,2);
    end
end
