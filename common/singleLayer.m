function S = singleLayer(rin,rout,mu)
%singleLayer(rin,rout,mu) determines the target-from-source matrix for
%2D Stokeslet sources centered at complex valued coordinates rin evaluated
%at targets rout, given viscosity mu
%
%code taken from Alex Barnett (BIE2D)

r = bsxfun(@minus, rout, rin.');          % C-# displacements mat
irr = 1./(conj(r).*r);                   % 1/r^2, used in all cases below
d1 = real(r); d2 = imag(r);              % worth storing I think
c = 1/(4*pi*mu);              % factor from Hsiao-Wendland book, Ladyzhenskaya

logir = -log(abs(r));  % log(1/r) diag block
A12 = d1.*d2.*irr;     % off diag vel block
S = c*[logir + d1.^2.*irr, A12;                         % u_x
     A12,                logir + d2.^2.*irr];         % u_y)
end