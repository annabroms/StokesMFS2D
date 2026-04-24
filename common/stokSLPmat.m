function S = stokSLPmat(rin,rout,mu)
%stokSLPmat(rin,rout,mu) determines the target-from-source matrix for
%2D Stokeslet sources centered at complex valued coordinates rin evaluated
%at targets rout, given viscosity mu
%
%code adapted from Alex Barnett (BIE2D)

rin = rin(:).';
rout = rout(:);

dx = real(rout) - real(rin);
dy = imag(rout) - imag(rin);
dx2 = dx.^2;
dy2 = dy.^2;
r2 = dx2 + dy2;
irr = 1./r2;

c = 1/(4*pi*mu); % factor from Hsiao-Wendland book, Ladyzhenskaya
logir = -0.5*log(r2); % log(1/r)
A12 = dx.*dy.*irr; % off diagonal velocity block
S = c*[logir + dx2.*irr, A12; ...
       A12,              logir + dy2.*irr];
end
