function L = lagpoly(n, alp)
%LAGPOLY   Chebfun representation of Laguerre polynomials.
%   L = LAGPOLY(N) returns the CHEBFUN corresponding to the Laguerre polynomials
%   L_N(x) on [0,inf]. N may be a vector of positive integers.
%
%   L = LAGPOLY(N, ALPA) is the same but for the generalized Laguerre
%   polynomials L^{(ALPHA)}_N(x).
%
%   Note, this is currently just a toy to play with the construction of Laguerre
%   polynomials using a combination of Chebfun's barycentric, mapping, and
%   'blowup' technologies.
%
% See also CHEBPOLY, LEGPOLY, JACPOLY, and HERMPOLY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    alp = 0;
end
dom = [0, inf];

if ( ~isscalar(alp) || ~isreal(alp) )
    error('CHEBFUN:chebfun:lapoly:alp', 'ALP must be a real-valued scalar');
end

x = chebfun(@(x) x, dom);                   % X on [0, inf]
L = chebfun(@(x) 1 + 0*x, dom);             % L_0(x)
L = [L, chebfun(@(x) 1 + alp - x, dom)];    % L_1(x)

for k = 2:max(n) % Recurrence relation
    Lk = ((2+(alp-1-x)/k).*L(:,k) - (1+(alp-1)/k)*L(:,k-1));
    L = [L, Lk]; %#ok<AGROW>
end

% Take only the ones we want:
L = L(:,n+1);

end