function p = ultrapoly(n, lam, dom)
%ULTRAPOLY   Ultraspherical polynomials.
%   P = ULTRAPOLY(N, LAM) computes a CHEBFUN of the ultraspherical polynomial of
%   degree N with parameters LAM, where the weight function is defined by w(x) =
%   (1-x^2)^(LAM-.5). N may be a vector of integers. LAM must be positive.
%
%   Normalization is chosen to be consistent with the formulas in [1, $18]. In
%   particular, feval(P, 1) = (2*LAM)_n/n! and feval(P, 0) = (-1)^n(LAM)_n/n!,
%   where ()_n is the Pochhammer notation for the rising factorial [1, (5.2.5)].
%
%   P = ULTRAPOLY(N, LAM, DOM) computes the ultraspherical polynomials as above,
%   but on the interval given by the domain DOM, which must be bounded.
%
%   P is computed via the standard recurrence relation for ultraspherical
%   polynomials.
%
%   References:
%    [1] F.W.J. Olver et al., editors. NIST Handbook of Mathematical Functions.
%    Cambridge University Press, New York, NY, 2010.
%
% See also LEGPOLY, CHEBPOLY, JACPOLY. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin < 2 )
    error('CHEBFUN:ultrapoly:inputs', 'ULTRAPOLY() requires at least 2 inputs.'); 
end
if ( nargin < 3 )
    dom = [-1, 1]; 
end
% Unbounded domains aren't supported/defined.
if ( any(isinf(dom)) )
    error('CHEBFUN:ultrapoly:infdomain', ...
        'Ultraspherical polynomials are not defined over an unbounded domain.');
end
if ( lam <= 0 )
    error('CHEBFUN:ultrapoly:infdomain', ...
        'Ultraspherical polynomials are not defined for LAMBA <= 0.');
end

if ( lam == .5 ) 
   % Equivalent to LEGPOLY: 
   p = legpoly(n, dom);
   return
elseif ( lam == 1 ) 
   % Equivalent to CHEBPOLY of 2nd kind:
   p = chebpoly(n, dom, 2);
   return
end

% Force a CHEBTECH basis.
defaultPref = chebfunpref();
pref = defaultPref;
tech = feval(pref.tech);
if ( ~isa(tech, 'chebtech') )
    pref.tech = @chebtech2;
end

% Construct the ultraspherical coefficients:
nn = 0:max(n);
scl = gamma(lam+.5)/gamma(2*lam)*exp(gammaln(2*lam+nn)-gammaln(lam+nn+.5));
c = diag(scl);
c = jac2cheb(c(:,n+1), lam - .5, lam - .5);

% Construct CHEBFUN from coeffs:
p = chebfun(c, dom([1, end]), pref, 'coeffs');   
p = restrict(p, dom);

% Adjust orientation:
if ( size(n, 1) > 1 )
   p = p.'; 
end

end