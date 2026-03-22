function [q,B] = grow_cluster(P,varargin)
%GROW_CLUSTER Grow a cluster of spheres/circles with a target nearest-neighbour gap.
%   q = GROW_CLUSTER(P,delta,dim,R,useKDTree,visualise,verbose)
%   builds a cluster of P centers in dim=2 or dim=3.
%   In 2D, q is returned as complex coordinates.
%
%   stats = GROW_CLUSTER('self_test',P,delta,dim,R,nTrials,seed,visualise,verbose)
%   runs brute-force vs KD-tree timing. nTrials is the number of repeated
%   timing runs used for averaging.
%
%   GROW_CLUSTER() with no inputs runs self_test with defaults. The kd-tree branch is 
%   slower unless P is large, and deactivated for small P.
%
%   Anna Broms 4/12/24, updated Mar 2026

if nargin == 0
    q = runSelfTest();
    B = [];
    return
end

if ischar(P) || (isstring(P) && isscalar(P))
    mode = lower(char(string(P)));
    if ~strcmp(mode,'self_test')
        error('Unknown mode "%s". Use numeric inputs or ''self_test''.',mode);
    end
    q = runSelfTest(varargin{:});
    B = [];
    return
end

if isempty(varargin)
    error('grow_cluster requires at least P and delta.');
end

if numel(varargin) > 6
    error('Too many inputs. Expected: P,delta,dim,R,useKDTree,visualise,verbose.');
end

delta = varargin{1};
if numel(varargin) >= 2 && ~isempty(varargin{2})
    dim = varargin{2};
else
    dim = 3;
end
if numel(varargin) >= 3 && ~isempty(varargin{3})
    R = varargin{3};
else
    R = 1;
end
if numel(varargin) >= 4 && ~isempty(varargin{4})
    useKDTree = varargin{4};
else
    useKDTree = [];
end
if numel(varargin) >= 5 && ~isempty(varargin{5})
    visualise = varargin{5};
else
    visualise = [];
end
if numel(varargin) >= 6 && ~isempty(varargin{6})
    verbose = varargin{6};
else
    verbose = [];
end

q = buildCluster(P,delta,dim,R,useKDTree,visualise,verbose);
B = [];

end


function q = buildCluster(P,delta,dim,R,useKDTree,visualise,verbose)
if ~(isscalar(P) && isnumeric(P) && isfinite(P) && P >= 1 && P == floor(P))
    error('P must be a positive integer scalar.');
end
if ~(isscalar(delta) && isnumeric(delta) && isfinite(delta) && delta > 0)
    error('delta must be a positive scalar.');
end
if ~(isscalar(dim) && (dim == 2 || dim == 3))
    error('dim must be 2 or 3.');
end
if ~(isscalar(R) && isnumeric(R) && isfinite(R) && R > 0)
    error('R must be a positive scalar.');
end

kdtreeAvailable = hasKDTreeSupport();
if isempty(useKDTree)
    kdAutoMinP = 5000;
    useKDTree = kdtreeAvailable && (P >= kdAutoMinP);
end
if isempty(visualise)
    visualise = (dim == 3);
end
if isempty(verbose)
    verbose = true;
end

useKDTree = asLogical(useKDTree,'useKDTree');
visualise = asLogical(visualise,'visualise');
verbose = asLogical(verbose,'verbose');

if useKDTree && ~kdtreeAvailable
    useKDTree = false;
    if verbose
        warning('grow_cluster:NoKDTree',...
            'KD-tree functions unavailable. Falling back to brute-force nearest-neighbour checks.');
    end
end

if dim == 3
    q = [0 0 0];
    hasLocalCenter = false;
else
    q = [0 0];
    hasLocalCenter = true;
end

if P >= 2
    n = randn(1,dim);
    n = n/norm(n);
    q(2,:) = R*(2+delta)*n;
end

nnTree = [];
treePts = zeros(0,dim);
bufferPts = zeros(0,dim);
if useKDTree
    % Delay first tree build to avoid overhead for small clusters, then
    % rebuild only in chunks while scanning a small append-buffer exactly.
    kdMinPoints = 400;
    kdBufferMax = 128;
    bufferPts = q;
end

maxit = 100;
tol = 1e-5;
deltaR = delta*R;

if verbose
    fprintf('== Generate particle cluster == \n');
end

for k = 1:P-2
    if verbose
        fprintf('Generate body %u/%u\n',k+2,P);
    end

    qExisting = q(1:k+1,:);
    useKDQuery = false;
    if useKDTree
        if isempty(nnTree) && size(qExisting,1) >= kdMinPoints
            treePts = qExisting;
            nnTree = createns(treePts,'NSMethod','kdtree','Distance','euclidean');
            bufferPts = zeros(0,dim);
        end
        useKDQuery = ~isempty(nnTree);
    end

    while true
        n = randn(1,dim);
        n = n/norm(n);

        if hasLocalCenter
            cn = randi(max(1,k-1));
            center = qExisting(cn,:);
        else
            center = [];
        end

        x1 = 2*R + deltaR + 0.1*R;
        x2 = 2*R + deltaR;
        f1 = minDist(qExisting,n,x1,center,R,nnTree,bufferPts,useKDQuery) - deltaR;
        f2 = minDist(qExisting,n,x2,center,R,nnTree,bufferPts,useKDQuery) - deltaR;
        itr = 0;

        while abs(f2)/delta > tol && itr < maxit
            denom = f2 - f1;
            if (denom == 0) || ~isfinite(denom)
                itr = maxit;
                break
            end
            xnew = x2 - f2*(x2-x1)/denom;
            if ~isfinite(xnew)
                itr = maxit;
                break
            end
            x1 = x2;
            f1 = f2;
            x2 = xnew;
            f2 = minDist(qExisting,n,x2,center,R,nnTree,bufferPts,useKDQuery) - deltaR;
            itr = itr + 1;
        end

        if itr < maxit
            break
        end
        if verbose
            disp('Cannot find config - try new direction');
        end
    end

    if hasLocalCenter
        q(k+2,:) = q(cn,:) + x2*n;
    else
        q(k+2,:) = x2*n;
    end

    if useKDTree
        if isempty(nnTree)
            bufferPts(end+1,:) = q(k+2,:);
        else
            bufferPts(end+1,:) = q(k+2,:);
            if size(bufferPts,1) >= kdBufferMax
                treePts = [treePts; bufferPts];
                nnTree = createns(treePts,'NSMethod','kdtree','Distance','euclidean');
                bufferPts = zeros(0,dim);
            end
        end
    end
end

if visualise
    visualizeCluster(q,dim,R,'Generated cluster');
end

if dim == 2
    q = q(:,1) + 1i*q(:,2);
end

if verbose
    disp('Cluster complete...');
end
end


function stats = runSelfTest(P,delta,dim,R,nTrials,seed,visualise,verbose)
if nargin < 1 || isempty(P)
    P = 10000;
end
if nargin < 2 || isempty(delta)
    delta = 0.001;
end
if nargin < 3 || isempty(dim)
    dim = 2;
end
if nargin < 4 || isempty(R)
    R = 1;
end
if nargin < 5 || isempty(nTrials)
    nTrials = 3;
end
if nargin < 6 || isempty(seed)
    seed = 1;
end
if nargin < 7 || isempty(visualise)
    visualise = false;
end
if nargin < 8 || isempty(verbose)
    verbose = true;
end
if nargin > 8
    error('Too many self_test inputs. Expected up to 8 optional values.');
end

visualise = asLogical(visualise,'visualise');
verbose = asLogical(verbose,'verbose');

if ~(isscalar(nTrials) && isnumeric(nTrials) && isfinite(nTrials) && nTrials >= 1 && nTrials == floor(nTrials))
    error('nTrials must be a positive integer.');
end

kdtreeAvailable = hasKDTreeSupport();
tBrute = zeros(nTrials,1);
tKD = nan(nTrials,1);
qBrute = [];
qKD = [];

if verbose
    fprintf('== grow_cluster self_test == \n');
    fprintf('P=%d, delta=%g, dim=%d, R=%g, nTrials=%d\n',P,delta,dim,R,nTrials);
end

for i = 1:nTrials
    rng(seed+i-1,'twister');
    t0 = tic;
    qBrute = buildCluster(P,delta,dim,R,false,false,false);
    tBrute(i) = toc(t0);

    if kdtreeAvailable
        rng(seed+i-1,'twister');
        t0 = tic;
        qKD = buildCluster(P,delta,dim,R,true,false,false);
        tKD(i) = toc(t0);
    end
end

if kdtreeAvailable
    qFinal = qKD;
    meanKD = mean(tKD);
    meanBrute = mean(tBrute);
    speedup = meanBrute/meanKD;
    maxDiff = max(abs(qBrute(:)-qKD(:)));
else
    qFinal = qBrute;
    meanKD = NaN;
    meanBrute = mean(tBrute);
    speedup = NaN;
    maxDiff = NaN;
end

stats = struct();
stats.P = P;
stats.delta = delta;
stats.dim = dim;
stats.R = R;
stats.nTrials = nTrials;
stats.seed = seed;
stats.kdTreeAvailable = kdtreeAvailable;
stats.bruteforceTimes = tBrute;
stats.kdTreeTimes = tKD;
stats.bruteforceMean = meanBrute;
stats.kdTreeMean = meanKD;
stats.speedupKDvsBruteforce = speedup;
stats.maxDiffSameSeed = maxDiff;
stats.minCenterDistance = minCenterDistance(qFinal,dim);

if verbose
    fprintf('Brute-force mean: %.6f s\n',meanBrute);
    if kdtreeAvailable
        fprintf('KD-tree mean:     %.6f s\n',meanKD);
        fprintf('Speedup:          %.2fx\n',speedup);
        fprintf('Max |q_brute-q_kd| (same seed): %.3e\n',maxDiff);
    else
        fprintf('KD-tree functions unavailable. KD timing skipped.\n');
    end
    fprintf('Minimum center distance in final cluster: %.6f\n',stats.minCenterDistance);
end

if visualise
    if kdtreeAvailable
        titleStr = sprintf('grow_cluster self_test (KD-tree, P=%d, dim=%d)',P,dim);
    else
        titleStr = sprintf('grow_cluster self_test (brute-force, P=%d, dim=%d)',P,dim);
    end
    visualizeCluster(qFinal,dim,R,titleStr);
end
end


function tf = hasKDTreeSupport()
tf = (exist('createns','file') == 2) && (exist('knnsearch','file') == 2);
end


function v = asLogical(v,name)
if ~(isscalar(v) && (islogical(v) || isnumeric(v)))
    error('%s must be a logical scalar.',name);
end
v = logical(v);
end


function dmin = minCenterDistance(q,dim)
if dim == 2
    pts = [real(q(:)), imag(q(:))];
else
    pts = q;
end

n = size(pts,1);
if n < 2
    dmin = inf;
    return
end

dmin = inf;
for i = 1:n-1
    dij = vecnorm(pts(i+1:end,:) - pts(i,:),2,2);
    dmin = min(dmin,min(dij));
end
end


function visualizeCluster(q,dim,R,titleStr)
if dim == 2
    if ~isreal(q)
        pts = [real(q(:)), imag(q(:))];
    else
        pts = q;
    end
    theta = linspace(0,2*pi,100);
    xc = R*cos(theta);
    yc = R*sin(theta);

    figure();
    hold on
    for k = 1:size(pts,1)
        plot(pts(k,1)+xc,pts(k,2)+yc,'LineWidth',0.8);
    end
    plot(pts(:,1),pts(:,2),'k.','MarkerSize',8);
    axis equal
    grid on
    xlabel('x');
    ylabel('y');
    title(titleStr);
else
    pts = q;
    [X,Y,Z] = sphere(15);

    figure();
    hold on
    for k = 1:size(pts,1)
        surf(pts(k,1)+R*X,pts(k,2)+R*Y,pts(k,3)+R*Z,...
            'EdgeColor','none','FaceAlpha',0.9);
    end
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(titleStr);
    colormap(copper);
    camlight;
    lighting gouraud;
end
end


function res = minDist(q,n,x,c,R,nnSearcher,bufferPts,useKDQuery)
if nargin < 6
    nnSearcher = [];
end
if nargin < 7
    bufferPts = zeros(0,size(q,2));
end
if nargin < 8
    useKDQuery = false;
end

if nargin < 4 || isempty(c)
    qnew = x*n;
else
    qnew = c + x*n;
end

if ~useKDQuery || isempty(nnSearcher)
    res = min(vecnorm(q-qnew,2,2)-2*R);
else
    [~,d] = knnsearch(nnSearcher,qnew,'K',1);
    if ~isempty(bufferPts)
        d = min(d,min(vecnorm(bufferPts-qnew,2,2)));
    end
    res = d - 2*R;
end
end
