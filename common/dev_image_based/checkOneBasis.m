function checkOneBasis(rbase_in_c,U,Y,L,Nc,a_c)
%Take U,Y as input, or recompute for an even coarser grid?

%% Evaluation points
%create grid to check basis in, at the boundaries of the particles in the
%pair
%parameterisation to visualise in at boundary
tself = linspace(0,2*pi,a_c*Nc+1); 
tself = tself(1:end-1)'; 
teval = linspace(0,2*pi,50*a_c*Nc+1); 
teval = teval(1:end-1)'; 

%% Tools needed to determine the pair-basis

%Initialize data on boundary
mu = zeros(2*a_c*Nc,1);

%Keep track of what coarse point on boundary that is nonzero
t = linspace(0,2*pi,a_c*Nc+1); 
t = t(1:end-1)'; 

%coarse self-interaction
Nself = stokSLPmat(rbase_in_c,cos(tself)+1i*sin(tself),1);
 

%% Check basis
% loop over x, 
for k = 1 :2*a_c*Nc 

    mu = cos(k*t); %different set of basis functions here
    mu = [mu; zeros(size(mu))];
    step1 = U{1}*mu; %Map to data at proxy points
    tau_mapped = Y{1}*step1;

    %If mobility, need to project
    if ~isempty(L)
        tau_mapped_non = tau_mapped(1:2*Nc);
        tau_i_x = tau_mapped(1:Nc);
        tau_i_y = tau_mapped(Nc+1:end);
        tau_mapped = [tau_i_x; tau_i_y]-L*[tau_i_x; tau_i_y];
    end

    %% Look first at the 1-body basis...
    % ...on the body itself. 
    uself = Nself*tau_mapped; 
    err(k) = norm(mu-uself);

%     figure(6);
%     clf;
%     plot(uself(1:end/2),'LineWidth',2)
%     hold on
%     plot(mu)


    
end


end

function removePatches(q,h,rad,maxval)
t = linspace(0,2*pi,200);
for k = 1:size(q,1)
    r_range = linspace(0,rad(k)+10*h,2);
    [R,T] = meshgrid(r_range,t);
    [X,Y] = pol2cart(T,R);
    patch(X(:)+q(k,1),Y(:)+q(k,2),maxval*ones(size(X(:))),'w','EdgeColor', 'w');
end

end