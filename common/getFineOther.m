function [rout_fine_other,tout] = getFineOther(a,N_f,varargin)
%Retrieve fine grid on circle p2 close to circle i

if numel(varargin) < 4
    error('getFineOther:BadInput','Expected at least 4 trailing inputs.');
end

rad = 1;
if iscell(varargin{1})
    refine = varargin{1};
    q = varargin{2};
    i = varargin{3}; %#ok<NASGU>
    p2 = varargin{4};
    if numel(varargin) >= 5 && ~isempty(varargin{5})
        rad = varargin{5};
    end
else
    rads = varargin{1};
    refine = varargin{2};
    q = varargin{3};
    i = varargin{4}; %#ok<NASGU>
    p2 = varargin{5};
    if isscalar(rads)
        rad = rads;
    else
        rad = rads(p2);
    end
end

nout = ceil(a*N_f); % number of uniformly distributed collocation points on fine grid.
t = linspace(0,2*pi,nout+1);
t = t(1:end-1)';
rout_base = rad*(cos(t)+1i*sin(t));

t2 = refine{p2,i}; %get parameterisation of the clustered points close to contact
fine_2 = q(p2)+rad*(cos(t2)+1i*sin(t2));

rout_fine_other = [q(p2)+rout_base; fine_2];

tout = [t; t2];


end
