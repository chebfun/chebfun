function varargout = gpr(x, y, varargin)
%GPR    Gaussian process regression
%
%   F = GPR(X, Y) returns a CHEBFUN F defined on [min(X),max(X)]
%   representing the posterior mean of a Gaussian process with prior mean 0
%   and squared exponential kernel
%               k(x,x') = SIGMA^2*exp(-1/(2*L^2)*(x-x')^2).
%   The default signal variance is SIGMA = max(abs(Y)). L is chosen such
%   that it maximizes the log marginal likelihood (see eq. (2.30) from [1]).
%   F matches Y at X.
% 
%   [F, FVAR] = GPR(X, Y) also returns a CHEBFUN representing an estimate
%   of the variance in the posterior.
%
%   [F, FVAR, SAMPLES] = GPR(X, Y, 'samples', N) also computes N
%   independent samples from the posterior distribution, returning them
%   as the N columns of the quasimatrix SAMPLES.
%
%   [...] = GPR(...,'domain', DOM) computes the results on the domain
%   DOM = [A, B].
%
%   [...] = GPR(...,'trig') uses a periodic version of the squared
%   exponential kernel (see eq. (4.31) from [1]), namely
%               k(x,x') = SIGMA^2*exp(-2/L^2*sin(pi*(x-x')/P)^2),
%   where P is the period length, corresponding to the size of the
%   approximation domain. If the domain is not specified, it is chosen
%   to be [min(X) max(X) +.1*(max(X)-min(X))] to account for the fact
%   that the Y values might not be identical at min(X) and max(X).
%
%   [...] = GPR(...,'sigma', SIGMA) specifies the signal variance of 
%   the kernel function.
%   [...] = GPR(...,'L', L) specifies the length scale parameter of the
%   kernel function.
%
%   [...] = GPR(...,'noise', sigmaY) specifies that the input is noisy
%   i.i.d. with noise distribution N(0,sigmaY^2). The kernel function
%   takes this into account by updating its values at the sample points. It
%   becomes
%           k'(x,x') = k(x,x') + sigmaY^2*delta_xx',
%   where delta_xx' is the Kronecker delta function taking one if x=x' and
%   zero otherwise (see eq. (2.20) from [1]).
%   
%
% Examples:
%
%      rng(1)
%      n = 10; x = -2 + 4*rand(n,1);
%      y = sin(exp(x));
%      [f,fvar,smpl] = gpr(x,y,'domain',[-2,2],'samples',3);
%      plot(f), hold on
%      plot(smpl,'color',[.8 .8 .8]), plot(x,y,'.k','markersize',14),
%      hold off
%      
%      % add noise to the inputs
%      y = y + .1*randn(n,1);
%      f2 = gpr(x,y,'domain',[-2,2],'noise',.1);
%      plot(f2), hold on
%      plot(x,y,'.k','markersize',14), hold off    
%
% Reference:
%
%   [1] C. E. Rasmussen & C. K. I. Williams, "Gaussian Processes
%   for Machine Learning", MIT Press, 2006
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

x = x(:); y = y(:);

scalingFactor = 1;
sigmaGiven = 0;

for k = 1:length(varargin)
    if ( strcmpi('sigma', varargin{k}) )
        sigmaGiven = 1;
    end
end

if ~isempty(y) && ~sigmaGiven
    scalingFactor = max(abs(y));
end

yn = y/scalingFactor;

opts = parseInputs(x, yn, varargin{:});
if ~opts.sigmaGiven
    opts.sigma = scalingFactor;
end

% Construct the kernel matrix corresponding to x. For the moment,
% we assume a Gaussian squared exponential kernel. (see for
% instance eq. (2.31) from [1])

if ~isempty(x)
    
    n = length(x);
    r = repmat(x,1,n) - repmat(x',n,1);
    if opts.trig
        K = opts.sigma^2*exp(-2/(opts.lenScale^2) * ...
                sin(pi/(opts.dom(end)-opts.dom(1))*r).^2) + ...
                opts.sigmaY^2*eye(n);
    else
        K = (opts.sigma^2)*exp(-1/(2*opts.lenScale^2)*r.^2) + ...
            opts.sigmaY^2*eye(n);
    end
    % compute the Cholesky decomposition of K
    if opts.sigmaY == 0
        L = chol(K+1e-15*scalingFactor^2*n*eye(n), 'lower');
    else
        L = chol(K+opts.sigmaY^2*eye(n), 'lower');
    end
    % coefficients of the radial basis function expansion of the mean
    alpha = L'\(L\y);

    % constuct a Chebfun approximation for the posterior distribution mean
    if opts.trig && ~opts.sigmaY
        f = chebfun(@(z) mean(alpha, x, z, opts), opts.dom, 'trig', ...
            'eps', 1e-12,'splitting','on');
    else
        f = chebfun(@(z) mean(alpha, x, z, opts), opts.dom, ...
            'eps', 1e-12,'splitting','on');
    end
                        
    % compute the predictive variance based on a large sample set
    sampleSize = min(20*n,2000);
    xSample = chebpts(sampleSize,opts.dom);
    % check if we need to add the noise term in the covariance matrix
    inX = ismember(xSample,x);
    % construct the input matrices
    rx = repmat(xSample,1,n) - repmat(x',sampleSize,1);
    rxs = repmat(xSample,1,sampleSize) - repmat(xSample',sampleSize,1);
    
    if opts.trig
        Ks = opts.sigma^2*exp(-2/(opts.lenScale^2) * ...
            sin(pi/(opts.dom(end)-opts.dom(1))*rx).^2) + ...
            opts.sigmaY^2*(rx == 0);
        
        Kss = opts.sigma^2*exp(-2/(opts.lenScale^2) * ...
            sin(pi/(opts.dom(end)-opts.dom(1)) * rxs).^2) + ...
            opts.sigmaY^2*diag(inX);
    else
        Ks = opts.sigma^2*exp(-1/(2*opts.lenScale^2) * rx.^2) + ...
            opts.sigmaY^2*(rx == 0);
            
        Kss = opts.sigma^2*exp(-1/(2*opts.lenScale^2) * rxs.^2) + ...
            opts.sigmaY^2*diag(inX);
    end

    v = L\(Ks');
                            
    fvar = spdiags(Kss - v'*v, 0);
    fvar = chebfun(fvar,opts.dom);
    
else % no data points given
    
    % we are assuming a zero mean on the prior
    f = chebfun(0,opts.dom);
    
    fvar = chebfun(opts.sigma^2,opts.dom);
end
fvar = simplify(fvar);

% Take samples from the posterior and construct Chebfun representations
% of them. For the moment, just sample at a large number of points and
% construct Chebfun representations.

if ( opts.samples > 0 )
    if ~isempty(x)
        Ls = chol(Kss - v'*v + 1e-12*scalingFactor^2*n*eye(sampleSize),...
            'lower');
        
        fSample = repmat(f(xSample), 1, opts.samples) + ...
                  Ls*randn(sampleSize, opts.samples);
    
        fSample = chebfun(fSample,opts.dom);
    else
        sampleSize = 1000;
        if opts.trig
            xSample = linspace(opts.dom(1),opts.dom(end),sampleSize)';
            rxs = repmat(xSample,1,sampleSize) - ...
                repmat(xSample',sampleSize,1);
            Kss = opts.sigmaf^2*exp(-2/(opts.lenScale^2) * ...
                sin(pi/diff(opts.dom) * rxs).^2);
            
            Ls = chol(Kss + 1e-12*scalingFactor^2*eye(sampleSize),'lower');
            
        else
            xSample = chebpts(sampleSize,opts.dom); 
            rxs = repmat(xSample,1,sampleSize) - ...
                repmat(xSample',sampleSize,1);
            Kss = (opts.sigmaf^2)*exp(-1/(2*opts.lenScale^2)* rxs.^2);
            
            Ls = chol(Kss + 1e-12*scalingFactor^2*eye(sampleSize),'lower');

        end
        
        fSample = repmat(f(xSample), 1, opts.samples) + ...
                        Ls*randn(sampleSize, opts.samples);
    
        if opts.trig
            fSample = chebfun(fSample,opts.dom,'trig');
        else
            fSample = chebfun(fSample,opts.dom);
        end
        
    end
    varargout = {f, fvar, fSample};
else
    varargout = {f, fvar};
end

end

function opts = parseInputs(x, y, varargin)

if length(x) ~= length(y)
    error('CHEBFUN:CHEBFUN:gpr:badInput', ...
             'The number of points and data values must be equal.');
end

opts.samples = 0;
opts.sigma = 0;
opts.sigmaGiven = 0;
opts.sigmaY = 0;
opts.lenScale = 0;
opts.dom = [];
opts.trig = 0;

for k = 1:length(varargin)
    if ( strcmpi('trig', varargin{k}) )
        opts.trig = k;
    end
end

if opts.trig
    varargin(opts.trig) = [];
end

for k = 1:2:length(varargin)
    if ( strcmpi('samples', varargin{k}) )
        opts.samples = varargin{k+1};
    elseif ( strcmpi('sigma', varargin{k}) )
        opts.sigmaGiven = 1;
        opts.sigma = varargin{k+1};
    elseif ( strcmpi('L', varargin{k}) )
        opts.lenScale = varargin{k+1};
    elseif ( strcmpi('domain', varargin{k}) )
        opts.dom = varargin{k+1};
    elseif ( strcmpi('noise', varargin{k}) )
        opts.sigmaY = varargin{k+1};
    else
        error('CHEBFUN:CHEBFUN:gpr:badInput', ...
            'Unrecognized sequence of input parameters.');
    end
end

if isempty(opts.dom) % domain not provided, default to [min(x) max(x)]
    if isempty(x)
        opts.dom = [-1 1];
    elseif length(x) == 1
        opts.dom = [x-1 x+1];
    elseif opts.trig
        diff = (max(x)-min(x))/10;
        opts.dom = [min(x) max(x)+diff];
    else
        opts.dom = [min(x) max(x)];
    end
end

if ~isempty(x) && opts.trig % if domain endpoints are among data points,
                            % check to see if periodicity is enforced
    [~,idMin] = min(x);
    [~,idMax] =  max(x);
    if opts.dom(1) == x(idMin) && opts.dom(end) == x(idMax)
        if y(idMin) ~= y(idMax)
        end
    end
end

if ~opts.lenScale % hyperparameters not specified
    n = length(x);
    if ~opts.sigmaGiven
        opts.sigma = 1;
    end
    % Construct a chebfun approximation of the log marginal likelihood
    % parametrized on the length scale. Use the length scale maximizing
    % this function.
    domSize = opts.dom(end)-opts.dom(1);
    if opts.trig
        searchDom = [1/(2*n) 10];
    else
        searchDom = [1/(2*pi*n)*domSize 10/pi*domSize];
    end
    
    % heuristic for reducing the optimization domain for the max log
    % marginal likelihood estimation
    fdom1 = logML(searchDom(1),x,y,opts);
    fdom2 = logML(searchDom(2),x,y,opts);
    while(fdom1 > fdom2 && searchDom(2)/searchDom(1) > 1+1e-4)
        newBound = searchDom(1)+(searchDom(2) - searchDom(1))/10;
        fdomnew = logML(newBound,x,y,opts);
        if (fdomnew > fdom1)
            break;
        else
            searchDom(2) = newBound;
            fdom2 = fdomnew;
        end
    end
    f = chebfun(@(z) logML(z,x,y,opts),searchDom, ...
        'eps',1e-6,'splitting','on');
    [~, opts.lenScale] = max(f);
end

end

% Computes the mean function estimate of the GP (using a Gaussian squared
% exponential kernel) at the points xEval
function fxEval = mean(alpha, x, xEval, opts)

n = length(x);
xEval = xEval(:);
m = length(xEval);
rx = repmat(xEval,1,n) - repmat(x',m,1);
if opts.trig
    Kss = opts.sigma^2*exp(-2/(opts.lenScale^2) * ...
        sin(pi/diff(opts.dom)*rx).^2) + opts.sigmaY*(rx == 0);
else
    Kss = opts.sigma^2*exp(-1/(2*opts.lenScale^2)*rx.^2) + ...
        opts.sigmaY*(rx == 0);
end

fxEval = Kss*alpha;

end

% Computes the log marginal likelihood estimate for a given array of
% hyperparameters (i.e., length scales)
function fxEval = logML(lenScale, x,y, opts)

fxEval = lenScale;
[r,c] = size(lenScale);
n = length(x);
rx = repmat(x,1,n) - repmat(x',n,1);
for i = 1:r
    for j = 1:c
        if opts.trig
            K = opts.sigma^2*exp(-2/(lenScale(i,j)^2) * ...
                    sin(pi/diff(opts.dom)*rx).^2);
        else
            K = opts.sigma^2*exp(-1/(2*lenScale(i,j)^2) * rx.^2);
        end
    
        % compute the Cholesky decomposition of K
        if opts.sigmaY ~= 0
            L = chol(K+opts.sigmaY^2*eye(n), 'lower');
        else
            L = chol(K+1e-15*n*opts.sigma^2*eye(n), 'lower');
        end
        alpha = L'\(L\y);
        % log marginal likelihood (see line 7 from Alg. 2.1 in [1])
        fxEval(i,j) = -.5*y'*alpha - trace(log(L)) - n/2*log(2*pi);
    end
end

end
