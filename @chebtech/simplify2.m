function f = simplify2(f, tol)
%SIMPLIFY  Remove small trailing Chebyshev coeffs of a happy CHEBTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  CHEBTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense. It does this by calling the routine STANDARDCHOP.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of EPS.  If
%  TOL is a row vector with as many columns as F, then TOL(k) will be used as
%  the simplification tolerance for column k of F.
%
% See also STANDARDCHOP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case.
if ( isempty(f) )
    return
end

% Do nothing to an unhappy CHEBTECH.
if ( ~f.ishappy )
    return
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
coeffs = f.coeffs;
[n, m] = size(coeffs);
coeffs = chebtech2.vals2coeffs(chebtech2.coeffs2vals(coeffs));

%%%% MODIFICATION FOR "EXTEND2"

v = chebtech2.coeffs2vals(coeffs);  % Chebyshev transform
ii = (1:n)';                        % indices of all points
di = [0;sign(randn(n-2,1));0];      % random perturbations of the indices
jj = ii + di;                       % indices perturbed
xx = chebpts(n);                    % relevant Cheb points
dx = diff(xx);                      % distances between Cheb pts
dd = 0*xx;                          % perturbation distances
ll = find(di<0);                    % indices that perturb to the left
dd(ll) = -dx(ll-1);                 % associated signed distance
rr = find(di>0);                    % indices that perturb to the right
dd(rr) = dx(rr);                    % associated distance
yy = v(jj)-v;                       % vertical jump with perturbation
ep = 0.5*2^(-52);                   % typical epsilon on [-1,1]
ww = 0*xx;                          % local weights
mm = find(di~=0);                   % indices that move
ww(mm) = ep./dd(mm);                % local weights
ww = ww.*rand(n,1);                 % randomized weight
v2 = v + ww.*yy;                    % perturbed data
coeffs = chebtech2.vals2coeffs(v2); % inverse Chebyshev transform

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the default tolerance if none was supplied.
if ( nargin < 2 )
    p = chebtech.techPref();
    tol = p.chebfuneps;
end

% Recast TOL as a row vector.k
if ( size(tol, 2) ~= m )
    tol = max(tol)*ones(1, m);
end

% Loop through columns to compute CUTOFF.
cutoff = 1;
for k = 1:m
    cutoff = max(cutoff, standardChop(coeffs(:,k), tol(k)));
end

% Take the minimum of CUTOFF and LENGTH(F). This is necessary when padding was
% required.
cutoff = min(cutoff, nold);

% Chop coefficients using CUTOFF.
f.coeffs = f.coeffs(1:cutoff,:);

end
