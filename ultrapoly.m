function p = ultrapoly(n, lam, dom)
%ULTRAPOLY   Ultraspherical polynomials.
%   P = ULTRAPOLY(N, LAM) computes a CHEBFUN of the ultraspherical polynomial of
%   degree N with parameters LAM, where the weight function is
%   defined by w(x) = (1-x^2)^(LAM-.5). N may be a vector of integers.
%
%   P = ULTRAPOLY(N, LAM, DOM) computes the ultraspherical polynomials as above,
%   but on the interval given by the domain DOM, which must be bounded.
%
%   P is computed via the standard recurrence relation for ultraspherical polynomials.
%
% See also LEGPOLY, CHEBPOLY, JACPOLY. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% Parse inputs:
if ( nargin < 2 )
    error('CHEBFUN:ultrapoly:inputs', 'ULTRAPOLY() requires at least 2 inputs.'); 
end
if ( nargin < 4 )
    dom = [-1, 1]; 
end
% Unbounded domains aren't supported/defined.
if ( any(isinf(dom)) )
    error('CHEBFUN:ultrapoly:infdomain', ...
        'Ultraspherical polynomials are not defined over an unbounded domain.');
end

if ( lam == .5 ) 
   % default to LEGPOLY: 
   p = legpoly(n, dom);
   return
end

if ( lam == 1 ) 
    % default to CHEBPOLY of 2nd kind: 
   p = chebpoly(n, dom, 2);
   return
end

%% Setup:

% Force a CHEBTECH basis.
defaultPref = chebfunpref();
pref = defaultPref;
tech = feval(pref.tech);
if ( ~isa(tech, 'chebtech') )
    pref.tech = @chebtech2;
end

% Useful values:
nMax = max(n);
nMax1 = nMax + 1;
domIn = dom;
dom = dom([1, end]);
x = chebpts(nMax1, 2);

%% Recurrence relation:

P = zeros(nMax1);
P(:,1) = 1;    
P(:,2) = 2*lam*x;   
for k = 1:nMax-1
    P(:,k+2) = 2*(k+lam)/(k+1)*x.*P(:,k+1) - (k+2*lam-1)/(k+1)*P(:,k);
end

%% Assemble output:

[ignored1, ignored2, cc] = unique(n);  % Extract required column indices.
cc = nMax1 + 1 - cc;                   % P is ordered low to high.
C = chebtech2.vals2coeffs(P(:,cc));    % Convert to coefficients
C = fliplr(C);                         % C is ordered low to high.

% Construct CHEBFUN from coeffs:
p = chebfun(C, dom, pref, 'coeffs');   

if ( numel(domIn) > 2)
    p = restrict(p, domIn);
end

% Adjust orientation:
if ( size(n, 1) > 1 )
   p = p.'; 
end

end