function Nimage = getImageKernels2D(rimage,nimage,rtest,mu,s)
%GETIMAGEKERNELS2D(rimage,nimage,rtest,mu,singvec) evaluates different
%fundamental solutions centered at rimage in target points rtest. mu is the
%viscosity and s a list of the singularities to include, with 
% s(1) stokeslets (1 or 0?)
% s(2) rotlets (1 or 0?)
% s(3) stresslets (1 or 0?)
% s(4) potential dipoles (1 or 0?)
%Potential dipole from Pozrikidis: Introduction to theoretical and
%computational fluid dynamics. 

% Stokeslets
if s(1)
    Nimage = singleLayer(rimage,rtest,mu);
else
    Nimage = [];
end

r = bsxfun(@minus, rtest, rimage.'); 
irr = 1./(conj(r).*r); %1/r^2
%irr = 1./abs(r).^2; %same thing
irrrr = irr.*irr;
d1 = real(r); d2 = imag(r); 

%Rotlets
if s(2)
    c = 1/(4*pi*mu);
    R = -c*[d2.*irr; -d1.*irr];
    Nimage = [Nimage R];
end

%Stresslets
if s(3)
    rdotny = real(conj(r).*nimage.');                % r.n_y
    rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
    T = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
         imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
    Nimage = [Nimage T];
end

%Stresslets
if s(5)
    nimage = ones(size(rimage)); %normal in x direction
    rdotny = real(conj(r).*nimage.');                % r.n_y
    rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
    T = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
         imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
    Nimage = [Nimage T];

    nimage = 1i*ones(size(rimage)); %normal in y direction
   % nimage = zeros(size(nimage));
    rdotny = real(conj(r).*nimage.');                % r.n_y
    rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
    T = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
         imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
    Nimage = [Nimage T];
end

% Potential diples
if s(4)
    c = 1/(2*pi);
    A12 = 2*d1.*d2.*irrrr;
    D = c*[-irr+2*d1.^2.*irrrr, A12; A12 -irr+2*d2.^2.*irrrr]; %-delta_ij/r^2+2\hat x_i\hat x_i/r^4
    Nimage = [Nimage D];
end





end