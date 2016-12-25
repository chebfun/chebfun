function f = simplify(f, tol)
%SIMPLIFY  Remove small trailing Fourier coeffs of a happy TRIGTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  TRIGTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense. It does this by calling the routine STANDARDCHOP.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of EPS.  If
%  TOL is a row vector with as many columns as F, then TOL(k) will be used as
%  the simplification tolerance for column k of F.
%
% See also STANDARDCHOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case.
if ( isempty(f) )
    return
end

% Do nothing to an unhappy TRIGTECH. 
if ( ~f.ishappy )
    return;
end

% STANDARDCHOP requires at least 17 coefficients to avoid outright rejection.
% STANDARDCHOP also employs a look ahead feature for detecting plateaus. For F
% with insufficient length the coefficients are padded using prolong. The
% following parameters are chosen explicitly to work with STANDARDCHOP. 
% See STANDARDCHOP for details.
nold = length(f);
N = max(17, round(nold*1.25 + 5));
f = prolong(f,N);

% After the coefficients of F have been padded with zeros an artificial plateau
% is created using the noisy output from the FFT. The slightly noisy plateau is
% required since STANDARDCHOP uses logarithms to detect plateaus and this has
% undesirable effects when the plateau is made up of all zeros.
coeffs = abs(f.coeffs(end:-1:1,:));
[n, m] = size(coeffs);
coeffs = trigtech.vals2coeffs(trigtech.coeffs2vals(coeffs));

% Use the default tolerance if none was supplied.
if ( nargin < 2 )
    p = trigtech.techPref();
    tol = p.chebfuneps;
end

% Recast TOL as a row vector.
if ( size(tol, 2) ~= m )
    tol = max(tol)*ones(1, m);
end

% In order to work with STANDARDCHOP, the coefficients of F are modified so that
% the entries corresponding to wave numbers k and -k appear sequentially in the
% new matrix of coefficients. These entries are also replaced by the sum of the
% absolute values of the k and -k coefficients.

% Need to handle odd/even cases separately.
isEven = mod(n, 2) == 0;
if ( isEven )
    coeffs = [coeffs(n,:) ; coeffs(n-1:-1:n/2+1,:) + coeffs(1:n/2-1,:) ; coeffs(n/2,:)];
else
    coeffs = [coeffs(n:-1:(n+1)/2+1,:) + coeffs(1:(n+1)/2-1,:) ; coeffs((n+1)/2,:)];
end
coeffs = flipud(coeffs);
coeffs = [coeffs(1,:) ; kron(coeffs(2:end,:),[1 ; 1])];

% Loop through columns to compute CUTOFF.
cutoff = 1;
for k = 1:m
    cutoff = max(cutoff, standardChop(coeffs(:,k), tol(k)));
end
cutoff = min(cutoff,nold);

% Divide CUTOFF by 2.
if ( mod(cutoff, 2) == 0 )
    cutoff = cutoff/2 + 1;
else
    cutoff = (cutoff - 1)/2 + 1;
end

% Now put the coefficients vector back together.
coeffs = f.coeffs;
if ( isEven )
    coeffs = [.5*coeffs(1,:) ; coeffs(2:n,:) ; .5*coeffs(1,:)];
    n = n + 1;
end

% Use CUTOFF to trim F.
mid = (n + 1)/2;
f.coeffs = coeffs(mid-cutoff+1:mid+cutoff-1,:);
f.values = trigtech.coeffs2vals(f.coeffs);

end
