function f = simplify(f, tol)
%SIMPLIFY  Remove small trailing Fourier coefficients of a happy TRIGTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  TRIGTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense: ||G - F|| < G.EPSLEVEL*G.VSCALE. It does this by removing
%  trailing coefficients of F that are relatively small; more precisely, those 
%  that are smaller in magnitude than the product of F.VSCALE and F.EPSLEVEL. 
%  G.EPSLEVEL is set to F.EPSLEVEL.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of 
%  F.EPSLEVEL as the relative threshold level for deciding whether a coefficient
%  is small enough to be removed. Here, G.EPSLEVEL is set to the maximum of 
%  F.EPSLEVEL and TOL.
%
% See also HAPPINESSCHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Do nothing to an unhappy TRIGTECH. ([TODO]: Is this the right thing to do?)
if ( ~f.ishappy )
    return;
end

% Grab the coefficients:
coeffs = abs(f.coeffs(end:-1:1,:));
[n,m] = size(coeffs);

% Need to handle odd/even cases separately.
isEven = ~mod(n, 2);
if isEven
    % In this case the negative cofficients have an additional term
    % corresponding to the cos(N/2*x) coefficient.
    coeffs = [coeffs(n,:);coeffs(n-1:-1:n/2+1,:)+coeffs(1:n/2-1,:);coeffs(n/2,:)];
else
    coeffs = [coeffs(n:-1:(n+1)/2+1,:)+coeffs(1:(n+1)/2-1,:);coeffs((n+1)/2,:)];
end
coeffs = flipud(coeffs);
if ( size(coeffs,1) < 17 )
    coeffs = [coeffs;zeros(17-size(coeffs,1),m)];
end

% Use the default tolerance if none was supplied:
if ( nargin < 2 )
    p = chebfunpref;
    tol = p.eps;
end
if length(tol) ~= m
    tol = max(tol)*ones(1,m);
end

% Loop through columns to compute cutoff
cutoff = 1;
for k = 1:m
    cutoff = max(cutoff,standardChop(coeffs(:,k),tol(k),1));
end

% Now put the coefficients vector back together.
coeffs = f.coeffs;
if ( ~mod(n,2) )
    coeffs = [.5*coeffs(n,:);coeffs(1:n-1,:);.5*coeffs(n,:)];
    n = n+1;
end
mid = (n+1)/2;
f.coeffs = coeffs(mid-cutoff+1:mid+cutoff-1,:);
f.values = trigtech.coeffs2vals(f.coeffs);

% Update epslevel:
f.epslevel = eps + 0*f.epslevel;

end
