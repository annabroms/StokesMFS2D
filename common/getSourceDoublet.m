function u = getSourceDoublet(sd, rimage,rcheck)
s = zeros(1,7);
s(7) = 1; 
Nim = getImageKernels2D(rimage,[],rcheck,1,s); 

u = Nim*sd;

end