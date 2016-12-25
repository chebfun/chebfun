function p = jacpoly(n, a, b, dom)
%JACPOLY   Jacobi polynomials.
%   P = JACPOLY(N, ALPHA, BETA) computes a CHEBFUN of the Jacobi polynomial of
%   degree N with parameters ALPHA and BETA, where the Jacobi weight function is
%   defined by w(x) = (1-x)^ALPHA*(1+x)^BETA. N may be a vector of integers.
%
%   Normalization is chosen to be consistent with the formulas in [1, $18]. In
%   particular, feval(P, 1) = (ALPHA+1)_n/n!, where ()_n is the Pochhammer
%   notation for the rising factorial [1, (5.2.5)].
%
%   P = JACPOLY(N, ALPHA, BETA, DOM) computes the Jacobi polynomials as above,
%   but on the interval given by the domain DOM, which must be bounded.
%
%   P is computed via the standard recurrence relation for Jacobi polynomials.
%
%   References:
%    [1] F.W.J. Olver et al., editors. NIST Handbook of Mathematical Functions.
%    Cambridge University Press, New York, NY, 2010.
%
% See also LEGPOLY, CHEBPOLY, ULTRAPOLY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Use QR to compute the values, as we do in LEGPOLY()?

% Parse inputs:
if ( nargin < 3 )
    error('CHEBFUN:jacpoly:inputs', 'JACPOLY() requires at least 3 inputs.'); 
end
if ( nargin < 4 )
    dom = [-1, 1]; 
end
% Unbounded domains aren't supported/defined.
if ( any(isinf(dom)) )
    error('CHEBFUN:jacpoly:infdomain', ...
        'Jacobi polynomials are not defined over an unbounded domain.');
end

% Force a CHEBTECH basis.
defaultPref = chebfunpref();
pref = defaultPref;
tech = feval(pref.tech);
if ( ~isa(tech, 'chebtech') )
    pref.tech = @chebtech2;
end

% Construct the Jacobi coefficients:
N = max(n) + 1;
c = eye(N);
c = jac2cheb(c(:,n+1), a, b);

% Construct a CHEBFUN from the coefficients:
p = chebfun(c, dom([1, end]), 'coeffs');
p = restrict(p, dom);

% Adjust orientation:
if ( size(n, 1) > 1 )
   p = p.'; 
end

end

