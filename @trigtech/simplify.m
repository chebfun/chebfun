function f = simplify(f, tol)
%SIMPLIFY  Remove small trailing Fourier coeffs of a happy TRIGTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  TRIGTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense. It does this by calling the routine STANDARDCHOP.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of
%  EPS. 
%
% See also STANDARDCHOP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case.
if ( isempty(f) )
    return
end

% Do nothing to an unhappy TRIGTECH. 
if ( ~f.ishappy )
    return;
end

% Grab the coefficients of F.
coeffs = abs(f.coeffs(end:-1:1,:));
[n, m] = size(coeffs);

% Use the default tolerance if none was supplied.
p = chebfunpref;
if ( nargin < 2 )
    tol = p.eps;
end

% Reshape TOL.
if ( size(tol, 2) ~= m )
    tol = max(max(tol),p.eps)*ones(1, m);
end

% In order to work with STANDARDCHOP, the coefficients of F are modified so that
% the entries corresponding to wave numbers k and -k appear sequentially in the
% new matrix of coefficients. These entries are also replaced by the sum of the
% absolute values of the k and -k coefficients.

% Need to handle odd/even cases separately.
isEven = mod(n, 2) == 0;
if isEven
    coeffs = [coeffs(n,:) ; coeffs(n-1:-1:n/2+1,:) + coeffs(1:n/2-1,:) ; coeffs(n/2,:)];
else
    coeffs = [coeffs(n:-1:(n+1)/2+1,:) + coeffs(1:(n+1)/2-1,:) ; coeffs((n+1)/2,:)];
end
coeffs = flipud(coeffs);
coeffs = [coeffs(1,:) ; kron(coeffs(2:end,:),[1 ; 1])];

% STANDARDCHOP requires at least 17 coefficients, so for F such that LENGTH(F) <
% 17, the coefficients are padded with entries between TOL^(7/6) and TOL.
% These parameters are chosen explicitly to work with STANDARDCHOP.
% See STANDARDCHOP for details.
N = max(17, round(n*1.25 + 5));
cfmins = min(abs(coeffs), [], 1);
cfmaxs = max(abs(coeffs), [], 1);
if ( n < N )
    coeffs = [coeffs ; ones(N - n, 1)* ...
              (max(tol.^(7/6), min(cfmins./cfmaxs,tol)).*cfmaxs)];
end

% Loop through columns to compute CUTOFF.
cutoff = 1;
for k = 1:m
    cutoff = max(cutoff, standardChop(coeffs(:,k), tol(k)));
end

% Divide CUTOFF by 2.
if ( mod(cutoff, 2) == 0 )
    cutoff = cutoff/2 + 1;
else
    cutoff = (cutoff - 1)/2 + 1;
end

% Now put the coefficients vector back together.
coeffs = f.coeffs;
if ( isEven )
    coeffs = [.5*coeffs(n,:) ; coeffs(1:n-1,:) ; .5*coeffs(n,:)];
    n = n + 1;
end

% Use CUTOFF to trim F.
mid = (n + 1)/2;
cutoff = min(cutoff, mid);
f.coeffs = coeffs(mid-cutoff+1:mid+cutoff-1,:);
f.values = trigtech.coeffs2vals(f.coeffs);

% Set F.EPSLEVEL to MATLAB EPS.
f.epslevel = eps + 0*f.epslevel;

end
