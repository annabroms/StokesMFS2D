function z = createPeanut(q1,q2,Np,debug)
%CREATEPEANUT creates discretisation of separation surface for two circular particles
%
% Syntax: z = createPeanut(q1,q2,Np,debug)
%
% Input: 
% q1    - Complex valued center coordinate of circle 1
% q2    - Complex valued center coordinate of circle 2
% Np    - Number of discrete points on peanut 
% debug - Boolean to draw peanut
%

delta = abs(q1-q2)-2;
theta = acos(real(q2-q1)/(2+delta));
if imag(q2)<imag(q1)
    theta = 2*pi-theta;
end
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];


alpha = acos((2+delta)/4);

% h1 = (2*pi-2*alpha)/Np/4;
% %t1 = alpha:h1:2*pi-alpha;
% t1 = linspace(alpha,2*pi-alpha,Np/4);
% %t2 = pi+alpha:h1:3*pi-alpha;
% t2 = linspace(pi+alpha,3*pi-alpha,Np/4);
% 
% h2 = (pi-2*alpha)/Np/4;
% %t3 = alpha+h2:h2:pi-alpha-h2;
% t3 = linspace(h1+alpha,pi-alpha-h1,Np/4);
% %t4 = pi+alpha+h2:h2:2*pi-alpha;
% t4 = linspace(pi+alpha+h1,2*pi-alpha-h1,Np/4);

%% The peanut consists of four segments
%h1 = (2*pi-2*alpha)/(3*Np)/8; %disc of particle curve 1
%t1 = alpha:h1:2*pi-alpha;
t1 = linspace(alpha,2*pi-alpha,3*Np/8);
%t2 = pi+alpha:h1:3*pi-alpha;
t2 = linspace(pi+alpha,3*pi-alpha,3*Np/8); %disc of particle curve 2
%t1 = pi; %debugging
%t2 = 0;
%warning('Only two points')

%For the remaining part of the peanut
%h2 = (pi-2*alpha)/Np/8; 
%t3 = alpha+h2:h2:pi-alpha-h2;
t3 = linspace(alpha,pi-alpha,Np/8+2);
t3 = t3(2:end-1);
%t4 = pi+alpha+h2:h2:2*pi-alpha;
t4 = linspace(pi+alpha,2*pi-alpha,Np/8+2);
t4 = t4(2:end-1); 


z1 = cos(t1)+1i*sin(t1);
z2 = 2+delta+cos(t2)+1i*sin(t2); 
z3 = 1+delta/2+cos(t3)+(-2*sin(alpha)+sin(t3))*1i;
z4 = 1+delta/2+cos(t4)+(2*sin(alpha)+sin(t4))*1i;

z = [z1 z2 z3 z4];
%z = [z1 z2];

%Now, rotate

zmap = R*[real(z); imag(z)];

zmap = zmap(1,:)'+1i*zmap(2,:)'; 

z = zmap-mean(zmap)+mean([q1,q2]);


if debug

    figure(101)
    subplot(1,2,1)
    plot(real(z1),imag(z1),'r.--');
    hold on
    plot(real(z2),imag(z2),'b+--');
    plot(real(z3),imag(z3),'mo--');
    plot(real(z4),imag(z4),'cs--');
    axis equal
    
    %subplot(1,2,2)
    figure(2)
    magenta = [0.8, 0.0, 0.8];
    plot(real(z),imag(z),'.','Color',magenta)
    axis equal
end


end


