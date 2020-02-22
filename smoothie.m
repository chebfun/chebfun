function f = smoothie(varargin)
%SMOOTHIE  Random function C-infinity but nowhere analytic
%   F = SMOOTHIE returns a CHEBFUN on [-1,1] corresponding to a 
%   function that is C-infinity but not analytic.  The function
%   is derived from a Fourier series with random coefficients of
%   root-exponentially decreasing amplitudes.
%   F is obtained by calling SMOOTHIE('trig') on an
%   interval 20% longer and restricting the result to [-1,1].
%
%   SMOOTHIE(DOM) returns a result with domain DOM = [A, B].
%   SMOOTHIE(N) returns a quasimatrix with N independent columns.
%   SMOOTHIE('trig') returns a periodic smoothie.
%   SMOOTHIE('complex') returns a complex smoothie.
%
% Examples
%
%   f = smoothie;
%   subplot(2,1,1), plot(f)
%   subplot(2,1,2), plotcoeffs(f)
%
%   plot(smoothie(2))
%
%   plot(smoothie('trig','complex')), axis equal
%
% See also RANDNFUN.

% Copyright 2020 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[n, dom, trig, cmplx] = parseInputs(varargin{:});
L = diff(dom);
m = round(ceil(2000*L)) + 1;

if cmplx   % complex case: call real case twice
    if trig
        f = ( smoothie(n, dom, 'trig') + 1i*smoothie(n, dom, 'trig') )/sqrt(2);
    else
        f = ( smoothie(n, dom) + 1i*smoothie(n, dom) )/sqrt(2);
    end
    return
end

if ~trig   % nonperiodic case: call periodic case and restrict
    dom2 = dom(1) + [0 1.2*diff(dom)];
    f = smoothie(n, dom2, 'trig');
    x = chebpts(round(2.5*m)+20, dom);  % this number is large enough...
    f = chebfun(f(x), dom);             % ...so this matches f{dom(1),dom(2)}
    f = simplify(f);   
    return
end

c = randn(m,n) + 1i*randn(m,n);         % random coeffs
c(1,:) = sqrt(2)*real(c(1,:));
c = (exp(-sqrt((1:m)'/L)).*c)/sqrt(L);  % root-expoential decay
c = [conj(c(end:-1:2,:)); c];           % symmetrize for real result

if cmplx
    error('CHEBFUN:smoothie:UnknownOption',...
        'Unknown input parameter.')
end

f = chebfun(c, dom, 'trig', 'coeffs');

end

function [n, dom, trig, cmplx] = parseInputs(varargin)

n = NaN;
dom = NaN;
trig = 0;
cmplx = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        if ( v(1) == 't' )
            trig = 1;
        elseif ( v(1) == 'c' )
            cmplx = 1;
        else
            error('CHEBFUN:smoothie','Unrecognized string input')
        end
    elseif ~isscalar(v)
        dom = v;
    else
        n = v;
    end
end

if isnan(n)
   n = 1;           % default number of columns
end
if isnan(dom)
   dom = [-1 1];    % default domain
end

end
