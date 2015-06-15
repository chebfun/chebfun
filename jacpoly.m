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

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Use QR to compute the values, as we do in LEGPOLY()?

%% Parse inputs:
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

% TODO: This could also be done using a weighed QR and JACPTS (see LEGPOLY).

%% Recurrence relation:

apb = a + b;
aa  = a * a;
bb  = b * b;
P = zeros(nMax1);
P(:,1) = 1;    
P(:,2) = 0.5*(2*(a + 1) + (apb + 2)*(x - 1));   
for k = 2:nMax
    k2 = 2*k;
    k2apb = k2 + apb;
    q1 =  k2*(k + apb)*(k2apb - 2);
    q2 = (k2apb - 1)*(aa - bb);
    q3 = (k2apb - 2)*(k2apb - 1)*k2apb;
    q4 =  2*(k + a - 1)*(k + b - 1)*k2apb;
    P(:,k+1) = ((q2 + q3*x).*P(:,k) - q4*P(:,k-1)) / q1;
end

%% Assemble output:
P = P(:,n+1);                    % Extract required columns
C = chebtech2.vals2coeffs(P);    % Convert to coefficients

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

