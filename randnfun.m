function f = randnfun(varargin)
%RANDNFUN   Random smooth function
%   F = RANDNFUN(LAMBDA) returns a smooth CHEBFUN on [-1,1] with maximum
%   frequency about 2pi/LAMBDA and standard normal distribution N(0,1)
%   at each point.  F can be regarded as a sample path of a Gaussian
%   process.  It is obtained by calling RANDNFUN(LAMBDA, 'trig') on an
%   interval 20% longer and restricting the result to [-1,1].
%
%   RANDNFUN(LAMBDA, DOM) returns a result with domain DOM = [A, B].
%
%   RANDNFUN(LAMBDA, N) returns a quasimatrix with N independent columns.
%
%   RANDNFUN(LAMBDA, 'norm') normalizes the output by dividing it by about
%   SQRT(2*LAMBDA), so white noise is approached in the limit LAMBDA -> 0.
%
%   RANDNFUN(LAMBDA, 'trig') returns a random periodic function.  This
%   is defined by a finite Fourier-Wiener series with independent normally
%   distributed coefficients of equal variance.
%
%   RANDNFUN(LAMBDA, 'complex') returns a complex random function.
%
%   RANDNFUN() uses the default value LAMBDA = 1.  Combinations such
%   as RANDNFUN(DOM) and RANDNFUN('norm', LAMBDA) are allowed so long as
%   N, if present, is preceded by an explicit specification of LAMBDA.
%
% Examples:
%
%   f = randnfun(0.1); std(f), plot(f)
%   plotcoeffs(f, '.'), xlim([0 200])
%
%   X = randnfun(.01,2); cov(X)
%
%   s = randi(100);
%   rng(s), f1 = randnfun(0.5,'norm',[0 10],3);
%   rng(s), f2 = randnfun(0.1,'norm',[0 10],3);
%   plot(cumsum(f1),'k',cumsum(f2),'r')
%
%   plot(cumsum(randnfun(.01,[0 5],'complex','norm'))), axis equal
%
% See also RANDNFUN2, RANDNFUNSPHERE, RANDNFUNDISK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[lambda, n, dom, normalize, trig, cmplx] = parseInputs(varargin{:});

if trig    % periodic case: finite Fourier-Wiener series.
           % See sec. 16.3 of Kahane, "Some Random Series of Functions".
           % Usually Fourier-Wiener series are written with sines and
           % cosines, but we use complex arithmetic because that's what
           % the trigtech constructor expects with 'coeffs'.
           % The reordered index vector ii is introduced so that
           % the random coefficients are chosen in the order
           % correponding to wave numbers 0, 1, -1, 2, -2, 3, ....

    L = diff(dom);
    m = round(L/lambda);
    c = randn(2*n, 2*m+1);                       % real numbers, var 1
    ii = [2*m+1:-2:1 2:2:2*m];
    c = c(:,ii)'; 
    c = (c(:,1:n) + 1i*c(:,n+(1:n)))/sqrt(2);    % complex coeffs, var 1
    if ~cmplx
        c = (c + flipud(conj(c)))/sqrt(2);       % real coeffs, var 1
    end
    if normalize
        c = c/sqrt(L/2);    % on [-1,1], coeffs from N(0,1) (finite F-W series)
    else
        c = c/sqrt(2*m+1);  % regardless of domain, function values from N(0,1)
    end
    f = chebfun(c, dom, 'trig', 'coeffs');

else       % nonperiodic case: call periodic case and restrict

    dom2 = dom(1) + [0 1.2*diff(dom)];
    m = round(diff(dom)/lambda);

    if lambda == inf    
        c = randn;
        if cmplx
            c = (c + 1i*randn)/sqrt(2)
        end
        f = chebfun(c, dom);        % random constant function
	return
    end

    if cmplx
        if normalize
            f = randnfun(lambda, n, dom2, 'norm', 'trig', 'complex');
        else
            f = randnfun(lambda, n, dom2, 'trig', 'complex');
        end
    else
        if normalize
            f = randnfun(lambda, n, dom2, 'norm', 'trig');
        else
            f = randnfun(lambda, n, dom2, 'trig');
        end
    end

           % restrict the result to the prescribed interval
 
    x = chebpts(5*m+5, dom);  % this number is large enough...
    f = chebfun(f(x), dom);   % ...so this is equiv. to f{dom(1),dom(2)}
    f = simplify(f, 1e-13);   % loosened tolerance gives clean Cheb series

end
end

function [lambda, n, dom, normalize, trig, cmplx] = parseInputs(varargin)

lambda = NaN;
n = NaN;
dom = NaN;
normalize = 0;
trig = 0;
cmplx = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        if ( v(1) == 'n' )
            normalize = 1;
        elseif ( v(1) == 't' )
            trig = 1;
        elseif ( v(1) == 'c' )
            cmplx = 1;
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

