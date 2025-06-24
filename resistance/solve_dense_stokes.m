function [x, reserr, coefnorm,solres] = solve_dense_stokes(rin, rout, rimage, nimage, s, b, regtol, scalecols,w)
%SOLVE_DENSE_STOKES  Solve a Stokes resistance BVP using direct pseudo-inverse for the
%entire system. 
%
%   [x, reserr, coefnorm,solres] = SOLVE_DENSE_STOKES(rin, rout, rimage, nimage, s, b, regtol, scalecols,w)
%
%   Constructs and solves a dense linear system for the 2D Stokes single-layer
%   potential representation, possibly including image contributions from
%   other fundamental solutions. 
%
%   INPUTS:
%     rin       - N array of complex valued interior source points (can be empty).
%     rout      - M array of complex valued boundary collocation points (targets).
%     rimage    - Array of clustered image source points (can be empty).
%     nimage    - Array of image normal directions (same size as rimage).
%     s         - 1 x 5 boolean array of to determine fundamental solutons at the 
%                 rimage nodes: [S R Tr, D,T], Stokeslets (S), Rotlets (R), Stresslet 
%                 with random "normal direction" (Tr), Potential dipoles (D), 
%                 Stresslets with directions ex, ey (T)
%     b         - 2M right-hand side vector (target velocities).
%     regtol    - Regularization parameter for truncated SVD
%     scalecols - Logical flag, whether to precondition from the right with
%                 column size of the system matrix.
%     w         - (optional) weights to be used in preconditioning from the
%                 left
%
%   OUTPUTS:
%     x         - Solution vector (density coefficients).
%     solres    - Residual Ax-b after regularization
%     reserr    - Residual error ‖A*x - b‖_∞ after regularization.
%     coefnorm  - Infinity norm of the solution vector ‖x‖_∞ 
%
%   NOTES:
%     - If column scaling is enabled, the solution is rescaled at the end.

    if nargin < 9
        w = ones(size(rout));
    end
                                         
    %% Prepare system matrix
    mu = 1; %viscosity
    if size(rin,1)
        Nio = singleLayer(rin,rout,mu);
    else
        Nio = [];
    end
    if size(rimage,1)
        Nimage = getImageKernels2D(rimage,nimage,rout,mu,s);
    else
        Nimage = [];
    end
    
    Ntot = [Nio Nimage];

    W = diag(repmat(sqrt(w),2,1));
    
    if scalecols
        col_norms = vecnorm(Ntot);
        D = diag(1 ./ col_norms);
        
        Ntot = Ntot*D;
        
    end
    NtotD = W*Ntot;
    b = W*b;

    %Possibly, do also preconditioning from the left
    [Y, U] = getPseudoFactors(NtotD, regtol, 0);
    x = Y*(U'*b);
    solres = NtotD*x-b;
    reserr = norm(solres,inf);
    coefnorm = norm(x,inf); %more reasonable to look at largest coeff 
    if scalecols
        x = D*x;
    end

end
    
