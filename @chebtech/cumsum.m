function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a CHEBTECH.
%   CUMSUM(F) is the indefinite integral of the CHEBTECH F with the constant of
%   integration chosen so that F(-1) = 0.
%
%   CUMSUM(F, 2) will take cumulative sum over the columns of F which is an
%   array-valued CHEBTECH.
%
% See also DIFF, SUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the CHEBTECH G of length n is represented as
%       \sum_{r=0}^{n-1} c_r T_r(x)
% its integral is represented with a CHEBTECH of length n+1 given by
%       \sum_{r=0}^{n} b_r T_r(x)
% where b_0 is determined from the constant of integration as
%       b_0 = \sum_{r=1}^{n} (-1)^(r+1) b_r;
% and other coefficients are given by
%       b_1 = c_0 - c_2/2,
%       b_r = (c_{r-1} - c_{r+1})/(2r) for r > 0,
% with c_{n+1} = c_{n+2} = 0.
%
% [Reference]: Pages 32-33 of Mason & Handscomb, "Chebyshev Polynomials".
% Chapman & Hall/CRC (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trivial case of an empty CHEBTECH:
if ( isempty(f) )
    return
end

% Sum with respect to the continuous variable by default:
if ( nargin < 2 )
    dim = 1;
end

if ( dim == 1 )
    % Take difference across 1st dimension:
    f = cumsumContinuousDim(f);
else
    % Take difference across 2nd dimension:
    f = cumsumFiniteDim(f);
end

end

function f = cumsumContinuousDim(f)
% CUMSUM over the continuous dimension.

% Initialise storage:
c = f.coeffs;                     % Obtain Chebyshev coefficients {c_r}

[n, m] = size(c);
c = [ zeros(2, m) ; c ];          % Pad with zeros
b = zeros(n-1, m);                % Initialize vector b = {b_r}

% Compute b_(n+1) ... b_2:
b(1:n-1,:) = (c(3:end-1,:) - c(1:end-3,:)) ./ repmat(2*(n:-1:2)', 1, m);
b(n,:) = c(end,:) - c(end-2,:)/2; % Compute b_1
v = ones(1, n);
v(end-1:-2:1) = -1;
b(n+1,:) = v*b;                   % Compute b_0 (satisfies f(-1) = 0)

% Recover coeffs:
f.coeffs = b;

% Update vscale: 
f.vscale = getvscl(f);

% Update epslevel:
f.epslevel = updateEpslevel(f);

% Simplify (as suggested in Chebfun ticket #128)
f = simplify(f);

% Ensure f(-1) = 0:
lval = get(f, 'lval');
f.coeffs(end,:) = f.coeffs(end,:) - lval;

end

function f = cumsumFiniteDim(f)
% CUMSUM over the finite dimension.

f.coeffs = cumsum(f.coeffs, 2);
newVscale = getvscl(f);
epslevelApprox = sum(f.epslevel.*f.vscale, 2)/sum(newVscale, 2); % TODO: Is this right?
f.epslevel = updateEpslevel(f, epslevelApprox);
f.vscale = newVscale;

end
