function Nimage = getImageKernels2D(rimage,nimage,rtest,mu,s)
%GETIMAGEKERNELS2D(rimage,nimage,rtest,mu,singvec) evaluates different
%fundamental solutions centered at rimage in target points rtest. mu is the
%viscosity and s a list of the singularities to include, with 
% s(1) stokeslets (1 or 0?)
% s(2) rotlets (1 or 0?)
% s(3) stresslets (1 or 0?) (one set direction)
% s(4) potential dipoles (1 or 0?)
% s(5) stresslets (1 or 0?) (two directions)
% s(6) Stokes doublets (1 or 0?) (two directions)
% s(7) Source doublets
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
if s(2) || s(6)
    c = 1/(4*pi*mu);
    R = -c*[d2.*irr; -d1.*irr];
    if s(2)
        Nimage = [Nimage R];
    end

end

%Stresslets
if s(3)
    rdotny = real(conj(r).*nimage.');                % r.n_y
    rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
    T = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
         imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
    Nimage = [Nimage T];
end

%Stresslets in two directions and/or Stokes doublets
if s(5) || s(6)
    nimage1 = ones(size(rimage)); %normal in x direction
    rdotny = real(conj(r).*nimage1.');                % r.n_y
    rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
    T1 = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
         imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
    
   

    nimage2 = 1i*ones(size(rimage)); %normal in y direction
    rdotny = real(conj(r).*nimage2.');                % r.n_y
    rdotnir4 = (1/pi) * rdotny.*irrrr;      % w/ overall prefac
    T2 = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
         imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
    if s(5)
        Nimage = [Nimage T1 T2];
    else
        SD1 = T1+R*[diag(imag(nimage1)) -diag(real(nimage1))]; %Stokes doublet
        SD2 = T2+R*[diag(imag(nimage2)) -diag(real(nimage2))]; %Stokes doublet (other direction)
        Nimage = [Nimage SD1 SD2];
    end
end

% Potential dipoles
if s(4)
    c = 1/(2*pi);
    A12 = 2*d1.*d2.*irrrr;
    D = c*[-irr+2*d1.^2.*irrrr, A12; A12 -irr+2*d2.^2.*irrrr]; %-delta_ij/r^2+2\hat x_i\hat x_i/r^4
    Nimage = [Nimage D];
end

% Source doublets (gradient of potential dipoles)
if s(7)
    c = 1/pi;
    
    % Common terms
    A = real(r); B = imag(r); % r = A + iB
    r4 = abs(r).^4;
    r6 = r4 .* abs(r).^2;

    % Direction 1: d = [1; 0]
    D1 = c * [...
        4*A.^3 ./ r6 - A ./ r4, 4*A.^2.*B ./ r6 - B ./ r4; 
        4*A.^2.*B ./ r6 - B ./ r4, 4*A.*B.^2 ./ r6 ];

    % Direction 2: d = [0; 1]
    D2 = c * [...
        4*A.^2.*B ./ r6, 4*A.*B.^2 ./ r6 - A ./ r4;
        4*A.*B.^2 ./ r6 - A ./ r4, 4*B.^3 ./ r6 - B ./ r4 ];

    % Append to Nimage
    Nimage = [Nimage D1 D2];
end






end