function [q,B] = grow_cluster(P,delta,dim,R)
%GROW_CLUSTER(P,delta,dim,R) grows a cluster of P spheres of radius R with each sphere
%delta*R away from at least one neighbour. dim can be set to 2 or 3 for 2D or 3D
% Anna Broms 4/12/24

if nargin<4
    R = 1;
elseif nargin<3
    R = 1;
    dim = 3;
    c = 0;
end

%set position of first sphere;
if dim == 3
    q = [0 0 0];
    c = 0;
else 
    c = 1; 
    q = [0 0];
end

%Pick random orientation for the second
n = randn(1,dim);
n = n/norm(n);


q(2,:) = R*(2+delta)*n; %store positions
maxit = 100;
%maxit = 50; 
tol = 1e-5;
visualise = 1; 

fprintf('== Generate particle cluster == \n')

for k = 1:P-2
    fprintf('Generate body %u/%u\n',k+2,P)
    
    while true

        %new direction vector 
        n = randn(1,dim);
        n = n/norm(n);
        if c
            %pick close neighbour
            cn = randi(max(1,k-1));
        end
        
    
        %want to solve for the distance from the origin
        x1 = 2*R+delta*R+0.1*R; %initial guess 
        %x1 = 2*delta*k*2;
        x2 = 2*R+delta*R;

        if c
            f = @(x) minDist(q,n,x,q(cn,:),R)-delta*R;
        else
    
            f = @(x) minDist(q,n,x,R)-delta*R;
        end
        itr = 0;
    
        while abs(f(x2))/delta>tol && itr < maxit
            xnew = x2 - f(x2)*(x2-x1)/(f(x2)-f(x1));
            x1 = x2;
            x2 = xnew;
            itr = itr +1;
            abs(f(x2));
        end

        if itr < maxit
            break;
        else
            disp('Cannot find config - try new direction')
        end
    end
    itrs(k) = itr;

    if c
        q(k+2,:) = q(cn,:)+x2*n;
    else
        q(k+2,:) = x2*n; 
    end
end

%visualise the particles
if visualise && (dim==3)
    [X,Y,Z] = sphere(15);
    
    figure()
    for k = 1:P
        surf(q(k,1)+X,q(k,2)+Y,q(k,3)+Z,'EdgeColor','flat');
        hold on
    end
    axis equal
    colormap(copper)
end

B = []; 
if dim == 2
    q = q(:,1)+1i*q(:,2); 
end

% %Debug
% % Check shortest distances
% for k = 1:P
%     d(k,:) = vecnorm(q(k,:)-q,2,2);
% 
%     %Smallest dist for each particle?
%     m = mink(d(k,:),2)-2*R;
%     mind(k) = m(2);
% end
% 
% 
% d = d(:);
% B = mink(d(d>0),5*P)-2*R; %Returns the np smallest distances
% 
% if min(B)<0
%     error('Overlap has occured')
% end
% 
disp('Cluster complete...')

end


function res = minDist(q,n,x,c,R)
%return shortest distance to any other particle
if nargin < 4
    qnew = x*n;
else
    qnew = c+x*n;
end

res = min(vecnorm(q-qnew,2,2)-2*R);

end
