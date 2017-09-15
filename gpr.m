function varargout = gpr(x, y, varargin)
%GPR        Gaussian Process regression
%
%   [MU, S2] = GPR(X, Y) returns a CHEBFUN on [min(X),max(X)] corresponding
%   to the posterior mean of a Gaussian Process with prior mean 0 and a
%   squared exponential kernel k(x,x') = SIGMAF^2*exp(-1/(2*L^2)*(x-x')^2),
%   with signal variance SIGMAF^2 = 1.21*max(abs(Y))^2 and length scale
%   L = 1/length(X). MU interpolates Y at X. S2 represents a chebfun
%   estimate of the variance in the posterior.
%
%   [MU, S2, SAMPLE] = GPR(X, Y, 'sample', N) also computes N samples from
%   the posterior distribution, returning them as N independent columns of
%   the quasimatrix SAMPLE.
%
%   [...] = GPR(...,'domain', DOM) computes the results on the domain
%   DOM = [A, B].
%
%   [...] = GPR(...,'trig') uses a periodic version of the squared
%   exponential kernel (see eq. (4.31) from [1]), namely
%               k(x,x') = SIGMAF^2*exp(-2/L^2*sin(pi*(x-x')/P)^2),
%   where P is the period length, corresponding to the size of the
%   approximation domain.
%
%   [...] = GPR(...,'hyperparam', [SIGMAF, L]) specifies the
%   hyperparameters of the kernel function.
%
% Example:
%
%   n = 10; x = -2 + 4*rand(n,1); x = sort(x);
%   y = sin(exp(x));
%   [mu,s2,smpl] = gpr(x,y,'domain',[-2,2],'sample',3,'hyperparams',[1,0.5]);
%   plot(repmat(mu,1,3)-smpl);
%   hold on, plot(x,zeros(n,1),'.k','markersize',14), hold off
%
% References:
%
%   [1] C. E. Rasmussen & C. K. I. Williams, "Gaussian Processes
%   for Machine Learning", MIT Press, 2006
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

x = x(:); y = y(:);
opts = parseInputs(x, y, varargin{:});

% Construct the kernel matrix corresponding to x. For the moment,
% we assume a Gaussian squared exponential kernel. (see for
% instance eq. (2.31) from [1])

if ~isempty(x)
    
    n = length(x);
    if opts.trig
        K = opts.sigmaf^2*exp(-2/(opts.lenScale^2) * ...
                sin(pi/(opts.dom(end)-opts.dom(1))*(repmat(x,1,n) - ...
                repmat(x',n,1))).^2);
    else
        K = (opts.sigmaf^2)*exp(-1/(2*opts.lenScale^2)*(repmat(x,1,n) - ...
                                repmat(x',n,1)).^2);
    end

    % compute the Cholesky decomposition of K
    L = chol(K+1e-15*n*eye(n), 'lower');
    % coefficients of the radial basis function expansion of the mean
    alpha = L'\(L\y);

    % constuct a Chebfun approximation for the posterior distribution mean
    if opts.trig
        mu = chebfun(@(z) mean(alpha, x, z, opts), opts.dom, 'trig', ...
            'eps', 1e-10);
    else
        mu = chebfun(@(z) mean(alpha, x, z, opts), opts.dom, 'eps', 1e-10);
    end
                        
    % compute the predictive variance based on a large sample set
    sampleSize = min(20*n,2000);
    xSample = chebpts(sampleSize,opts.dom);
    
    if opts.trig
        Ks = opts.sigmaf^2*exp(-2/(opts.lenScale^2) * ...
            sin(pi/(opts.dom(end)-opts.dom(1))*(repmat(xSample,1,n) - ...
            repmat(x',sampleSize,1))).^2);
        
        Kss = opts.sigmaf^2*exp(-2/(opts.lenScale^2) * ...
            sin(pi/(opts.dom(end)-opts.dom(1)) * ...
            (repmat(xSample,1,sampleSize) - ...
            repmat(xSample',sampleSize,1))).^2);
    else
        Ks = opts.sigmaf^2*exp(-1/(2*opts.lenScale^2) * ...
                (repmat(xSample,1,n)-repmat(x',sampleSize,1)).^2);
            
        Kss = opts.sigmaf^2*exp(-1/(2*opts.lenScale^2) * ...
            (repmat(xSample,1, sampleSize) - ...
            repmat(xSample',sampleSize,1)).^2);
    end

    v = L\(Ks');
                            
    s2 = spdiags(Kss - v'*v, 0);
    s2 = chebfun(s2,opts.dom);
    
else % no data points given
    
    % we are assuming a zero mean on the prior
    mu = chebfun(0,opts.dom);
    
    s2 = chebfun(opts.sigmaf^2,opts.dom);
end

% Take samples from the posterior and construct Chebfun representations
% of them. For the moment, just sample at a large number of points and
% construct Chebfun representations.
if ( opts.samples > 0 )
    if ~isempty(x)
        Ls = chol(Kss - v'*v + 1e-12*n*eye(sampleSize),'lower');
        
        fSample = repmat(mu(xSample), 1, opts.samples) + ...
                  Ls*randn(sampleSize, opts.samples);
    
        fSample = chebfun(fSample,opts.dom);
    else
        sampleSize = 1000;
        if opts.trig
            xSample = linspace(opts.dom(1),opts.dom(end),sampleSize)';
            Kss = opts.sigmaf^2*exp(-2/(opts.lenScale^2) * ...
                sin(pi/(opts.dom(end)-opts.dom(1)) * ...
                (repmat(xSample,1,sampleSize) - ...
                repmat(xSample',sampleSize,1))).^2);
            
            Ls = chol(Kss + 1e-12*eye(sampleSize),'lower');
            
        else
            xSample = chebpts(sampleSize,opts.dom);          
            Kss = (opts.sigmaf^2)*exp(-1/(2*opts.lenScale^2)* ...
                (repmat(xSample,1, sampleSize) - ...
                repmat(xSample',sampleSize,1)).^2);
            
            Ls = chol(Kss + 1e-12*eye(sampleSize),'lower');
        end
        
        fSample = repmat(mu(xSample), 1, opts.samples) + ...
                        Ls*randn(sampleSize, opts.samples);
    
        if opts.trig
            fSample = chefbun(fSample,opts.dom,'trig');
        else
            fSample = chebfun(fSample,opts.dom);
        end
        
    end
    varargout = {mu, s2, fSample};
else
    varargout = {mu, s2};
end

end

function opts = parseInputs(x, y, varargin)

if length(x) ~= length(y)
    error('CHEBFUN:CHEBFUN:gpr:badInput', ...
             'The number of points and data values must be equal.');
end

opts.samples = 0;
opts.sigmaf = 0;
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
    if ( strcmpi('sample', varargin{k}) )
        opts.samples = varargin{k+1};
    elseif ( strcmpi('hyperparams', varargin{k}) )
        hyperparams = varargin{k+1};
        opts.sigmaf = hyperparams(1);
        opts.lenScale = hyperparams(end);
    elseif ( strcmpi('domain', varargin{k}) )
        opts.dom = varargin{k+1};
    else
        error('CHEBFUN:CHEBFUN:gpr:badInput', ...
            'Unrecognized sequence of input parameters.');
    end
end

if isempty(opts.dom) % domain not provided, default to [min(x) max(x)]
    opts.dom = [min(x) max(x)];
end

if opts.trig % if domain endpoints are among data points, check to see if
             % periodicity is enforced
             % TODO: allow some tolerences?
    [~,idMin] = min(x);
    [~,idMax] =  max(x);
    if opts.dom(1) == x(idMin) && opts.dom(end) == x(idMax)
        if y(idMin) ~= y(idMax)
        end
    end
end

if ~opts.sigmaf && ~opts.lenScale % hyperparameters not specified; for the
                                  % moment, just give some heuristic values
                                  % based on the input data
    n = length(x);
    opts.sigmaf = 1.1*sqrt(max(abs(y)));
    opts.lenScale = 1/n;
end

end

% Computes the mean function estimate of the GP (using a Gaussian squared
% exponential kernel)
function fxEval = mean(alpha, x, xEval, opts)

n = length(x);
xEval = xEval(:);
m = length(xEval);
if opts.trig
    Kss = opts.sigmaf^2*exp(-2/(opts.lenScale^2) * ...
        sin(pi/(opts.dom(end)-opts.dom(1))*(repmat(xEval,1,n) - ...
        repmat(x',m,1))).^2);
else
    Kss = opts.sigmaf^2*exp(-1/(2*opts.lenScale^2)*(repmat(xEval,1,n) - ...
                                            repmat(x',m,1)).^2);
end

fxEval = Kss*alpha;

end