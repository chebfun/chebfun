function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a CHEBTECH.
%   CUMSUM(F) is the indefinite integral of the CHEBTECH F with the constant of
%   integration chosen so that F(-1) = 0.
%
%   CUMSUM(F, 2) will take cumulative sum over the columns of F which is an
%   array-valued CHEBTECH.
%
% See also DIFF, SUM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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
    % cumsum across 1st dimension:
    f = cumsumContinuousDim(f);
else
    % cumsum across 2nd dimension:
    f.coeffs = cumsum(f.coeffs, 2);
end

end

function f = cumsumContinuousDim(f)
% CUMSUM over the continuous dimension.

% Initialise storage:
c = f.coeffs;                      % Obtain Chebyshev coefficients {c_r}

[n, m] = size(c);
c = [ c ; zeros(2, m) ;];          % Pad with zeros
b = zeros(n+1, m);                 % Initialize vector b = {b_r}

% Compute b_(2) ... b_(n+1):
b(3:n+1,:) = (c(2:n,:) - c(4:n+2,:)) ./ repmat(2*(2:n).', 1, m);
b(2,:) = c(1,:) - c(3,:)/2;        % Compute b_1
v = ones(1, n);
v(2:2:end) = -1;
b(1,:) = v*b(2:end,:);             % Compute b_0 (satisfies f(-1) = 0)

% Recover coeffs:
f.coeffs = b;

% Simplify (as suggested in Chebfun ticket #128)
f = simplify(f);

% Ensure f(-1) = 0:
lval = get(f, 'lval');
f.coeffs(1,:) = f.coeffs(1,:) - lval;

end
