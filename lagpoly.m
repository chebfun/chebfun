function L = lagpoly(n)
%LAGPOLY   Chebfun representation of Laguerre polynomials.
%   L = LAGPOLY(N) returns the CHEBFUN corresponding to the Laguerre polynomials
%   L_N(x) on [0,inf]. N may be a vector of positive integers.
%
%   Note, this is currently just a toy to play with the construction of Laguerre
%   polynomials using a combination of Chebfun's barycentric, mapping, and
%   'blowup' technologies.
%
% See also CHEBPOLY, LEGPOLY, JACPOLY, and HERMPOLY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

x = chebfun(@(x) x, [0, inf]);          % X on [0, inf]
L = chebfun(@(x) 1 + 0*x, [0, inf]);    % L_0(x)
L = [L, chebfun(@(x) 1 - x, [0, inf])]; % L_1(x)

for k = 2:max(n) % Recurrence relation
    Lk = ((2*k-1-x).*L(:,k) - (k-1)*L(:,k-1))/k;
    L = [L, Lk]; %#ok<AGROW>
end

% Take only the ones we want:
L = L(:,n+1);

end
