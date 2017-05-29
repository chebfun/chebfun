function f = randnfun(varargin)
%RANDNFUN   Random smooth function
%   F = RANDNFUN(LAMBDA) returns a smooth CHEBFUN on [-1,1] with maximum
%   frequency about 2pi/LAMBDA and standard normal distribution N(0,1)
%   at each point.  F can be regarded as one sample path of a Gaussian
%   process.  It is obtained by calling RANDNFUN(LAMBDA, 'trig') on an
%   interval about 20% longer and restricting the result to [-1,1].
%
%   RANDNFUN(LAMBDA, DOM) returns a result with domain DOM = [A, B].
%
%   RANDNFUN(LAMBDA, N) returns a quasimatrix with N independent columns.
%
%   RANDNFUN(LAMBDA, 'norm') normalizes the output by dividing it by
%   SQRT(LAMBDA), so white noise is approached in the limit LAMBDA -> 0.
%
%   RANDNFUN(LAMBDA, 'trig') returns a random periodic function.  This
%   is defined by a finite Fourier series with independent normally
%   distributed coefficients of equal variance.
%
%   RANDNFUN() uses the default value LAMBDA = 1.  Combinations such
%   as RANDNFUN(DOM), RANDNFUN('norm', LAMBDA) are allowed so long as
%   LAMBDA, if present, precedes N, if present.
%
% Examples:
%
%   f = randnfun(0.1); std(f), plot(f)
%   plotcoeffs(f, '.'), xlim([0 200])
%
%   X = randnfun(.01,2); cov(X)
%
%   f = randnfun(0.1,'norm',[0 10]); plot(cumsum(f))
%
% See also RANDNFUN2, RANDNFUNSPHERE, RANDNFUNDISK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[lambda, n, dom, normalize, trig] = parseInputs(varargin{:});

if trig    % periodic case: random Fourier series

    m = round(diff(dom)/lambda);
    c = randn(2*m+1, n) + 1i*randn(2*m+1, n);
    c = (c + flipud(conj(c)))/2;
    f = chebfun(c/sqrt(2*m+1), dom, 'trig', 'coeffs');
    if normalize
        f = f/sqrt(lambda);
    end

else       % nonperiodic case: call periodic case and restrict

    if lambda == inf    
        f = chebfun(randn, dom);        % random constant function
        if normalize
            f = f/sqrt(lambda);         % zero function
        end
        return
    end

    m = round(1.2*diff(dom)/lambda+2);
    dom2 = dom(1) + [0 m*lambda];
    if normalize
        f = randnfun(lambda, n, dom2, 'norm', 'trig');
    else
        f = randnfun(lambda, n, dom2, 'trig');
    end

           % restrict the result to the prescribed interval

    x = chebpts(5*m, dom);    % 5*m is large enough...
    f = chebfun(f(x), dom);   % ...so this is equiv. to f{dom(1),dom(2)}
    f = simplify(f, 1e-13);   % loosened tolerance gives clean Cheb series

end
end

function [lambda, n, dom, normalize, trig] = parseInputs(varargin)

lambda = NaN;
n = NaN;
dom = NaN;
normalize = 0;
trig = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        if ( v(1) == 'n' )
            normalize = 1;
        elseif ( v(1) == 't' )
            trig = 1;
        else
            error('CHEBFUN:randnfun','Unrecognized string input')
        end
    elseif ~isscalar(v)
        dom = v;
    elseif isnan(lambda)
        lambda = v;
    else
        n = v;
    end
end

if isnan(lambda)
   lambda = 1;      % default space scale
end
if isnan(n)
   n = 1;           % default number of columns
end
if isnan(dom)
   dom = [-1 1];    % default domain
end

end

