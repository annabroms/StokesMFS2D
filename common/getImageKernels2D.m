function Nimage = getImageKernels2D(rimage,nimage,rtest,mu,singvec)
%GETIMAGEKERNELS2D(rimage,nimage,rtest,mu,singvec) evaluates different
%fundamental solutions centered at rimage in target points rtest. mu is the
%viscosity and singvec a list of the singularities to include, with 
% s(1) stokeslets
% s(2) rotlets
% s(3) stresslets
% s(4) potential dipoles

%Potential dipole from Pozrikidis: Introduction to theoretical and
%computational fluid dynamics. 

% r = bsxfun(@minus, rtest, rimage.'); 
% irr = 1./(conj(r).*r); %1/r^2
% irrrr = irr.*irr;
% d1 = real(r); d2 = imag(r); 
% c = 1/(2*pi);
% A12 = 2*d1.*d2.*irrrr;
% D = c*[-irr+2*d1.^2.*irrrr, A12; A12 -irr+2*d2.^2.*irrrr]; %-delta_ij/r^2+2\hat x_i\hat x_i/r^4
% 
% %Rotlet
% %c = 1/(4*pi*mu);
% %R = -c*(d2 - 1i*d1)./irr;
% %R = -c*[d2./irr; -d1./irr];
% 
% %R = -c*[d2./irr -d1./irr; -d1./irr d2./irr]; %is this reasonable??
% %Stresslet
% rdotny = real(conj(r).*nimage.');                % r.n_y
% rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
% S = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
%      imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
% 


%Potential dipole from Pozrikidis: Introduction to theoretical and
%computational fluid dynamics. 

r = bsxfun(@minus, rtest, rimage.'); 
irr = 1./(conj(r).*r); %1/r^2
%irr_b = 1./abs(r).^2; %same thing
irrrr = irr.*irr;
d1 = real(r); d2 = imag(r); 
if singvec(4)
    c = 1/(2*pi);
    A12 = 2*d1.*d2.*irrrr;
    D = c*[-irr+2*d1.^2.*irrrr, A12; A12 -irr+2*d2.^2.*irrrr]; %-delta_ij/r^2+2\hat x_i\hat x_i/r^4
end

%Rotlet
if singvec(2)
    c = 1/(4*pi*mu);
    %R = -c*(d2 - 1i*d1)./irr;
    %R = [-c*d2./irr; c*d1./irr];
    R = -c*[d2.*irr; -d1.*irr];
end


%Stresslet
if singvec(3)
    rdotny = real(conj(r).*nimage.');                % r.n_y
    rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
    S = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
         imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
end

if singvec(1)
    Nimage = singleLayer(rimage,rtest,mu);
else
    Nimage = [];
end


if singvec(2)
    Nimage = [Nimage R];
end


if singvec(3)
    Nimage = [Nimage S];
end

if singvec(4)
    Nimage = [Nimage D];
end


end