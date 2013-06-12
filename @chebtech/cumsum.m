function f = cumsum(f, m, pref)
%CUMSUM   Indefinite integral of a CHEBTECH.
%   CUMSUM(F) is the indefinite integral of the CHEBTECH F with the constant of
%   integration chosen so that F(-1) = 0. 
%
%   CUMSUM(F, M) will compute the Mth definite integral with the constant of
%   integration chosen so that each intermediary integral evaluates to 0 at -1.
%   Thus, CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, PREF) or CUMSUM(F, M,  PREF) uses options from the preference
%   structure PREF when building the output CHEBTECH.
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

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

% Parse inputs:
if ( nargin == 1 )
    m = 1;
    pref = chebtech.pref();
elseif ( nargin < 3 )
    if ( isstruct(m) )
        pref = m;
        m = 1;
    else
        pref = chebtech.pref();
    end
end

% Initialise storage:
c = f.coeffs;                         % Obtain Chebyshev coefficients c = {c_r}

% Loop for higher-order integrals:
for k = 1:m

    [n, m] = size(c);
    c = [ zeros(2,m) ; c ];           %#ok<AGROW> % Pad with zeros
    b = zeros(n-1, m);                % Initialize vector b = {b_r}

    % Compute b_(n+1) ... b_2:
    b(1:n-1,:) = (c(3:end-1,:) - c(1:end-3,:)) ./ repmat(2*(n:-1:2)', 1, m);
    b(n,:) = c(end,:) - c(end-2,:)/2; % Compute b_1
    v = ones(1,n);
    v(end-1:-2:1) = -1;
    b(n+1,:) = v*b;                   % Compute b_0 (satisfies f(-1) = 0)
    
    % Copy coefficients back into c for loop:
    c = b;
end

% Recover values and attach to output:
f.values = f.chebpolyval(c);
f.coeffs = c;

% Update vscale:
f.vscale = max(f.vscale, max(abs(f.values), [], 1));
% [TODO]: Update epslevel?

% Simplify (as suggested in Chebfun ticket #128)
f = simplify(f);

% Ensure f(-1) = 0:
lval = get(f, 'lval');
f.coeffs(end,:) = f.coeffs(end,:) - lval;
f.values = bsxfun(@minus, f.values, lval);

end
