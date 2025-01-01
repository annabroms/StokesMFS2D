function qnew = createDumbells(P,delta)
%Create geometry with pairs of circles separated by delta. No triplets close. 
% P assumed to be even

if (P/2)>floor(P/2)
    error('P must be even')
end

[q,~] = grow_cluster(P/2,2+2*delta,2); %pick some separation between the dumbells here
%[q,B] = grow_cluster(P/2,2+1,2); %Testing
qnew = zeros(P,1); 
disp('Dumbell positions generated')
for k = 1:P/2
    %create random rotation 
    theta = 2*pi*rand(1);
    qtemp = (2+delta)*exp(1i*theta);
    q1 = -qtemp/2;
    q2 = qtemp/2;
    
    %also need displacement from center of big circle     
    qnew(2*(k-1)+1) = q(k)+q1;
    qnew(2*k) = q(k)+q2;

end



end