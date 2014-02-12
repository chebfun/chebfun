function f = cumsum(f, m, dim)
%CUMSUM   Indefinite integral of a CHEBTECH.
%   CUMSUM(F) is the indefinite integral of the CHEBTECH F with the constant of
%   integration chosen so that F(-1) = 0. 
%
%   CUMSUM(F, M) will compute the Mth definite integral with the constant of
%   integration chosen so that each intermediary integral evaluates to 0 at -1.
%   Thus, CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, M, 2) will take the Mth cumulative sum over the columns F an
%   array-valued CHEBTECH.
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

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

if ( nargin < 2 || isempty(m) )
    % Order of intergration not passed in. Assume 1 by default:
    m = 1; 
elseif ( m == 0 )
    % Nothing to do here!
    return
end    

% Sum with respect to the continuous variable by default:
if ( nargin < 3 )
    dim = 1;
end

if ( dim == 1 )
    % Take difference across 1st dimension:
    f = cumsumContinuousDim(f, m);
else
    % Take difference across 2nd dimension:
    f = cumsumFiniteDim(f, m);
end

end

function f = cumsumContinuousDim(f, m)
% CUMSUM over the continuous dimension.

    % Initialise storage:
    c = f.coeffs;                         % Obtain Chebyshev coefficients {c_r}

    % Loop for higher-order integrals:
    for k = 1:m

        [n, m] = size(c);
        c = [ zeros(2, m) ; c ];          %#ok<AGROW> % Pad with zeros
        b = zeros(n-1, m);                % Initialize vector b = {b_r}

        % Compute b_(n+1) ... b_2:
        b(1:n-1,:) = (c(3:end-1,:) - c(1:end-3,:)) ./ repmat(2*(n:-1:2)', 1, m);
        b(n,:) = c(end,:) - c(end-2,:)/2; % Compute b_1
        v = ones(1, n);
        v(end-1:-2:1) = -1;
        b(n+1,:) = v*b;                   % Compute b_0 (satisfies f(-1) = 0)

        % Copy coefficients back into c for loop:
        c = b;
    end

    % Recover values and attach to output:
    f.values = f.coeffs2vals(c);
    f.coeffs = c;

    % Update vscale: [TODO]: Update epslevel?
    f.vscale = max(abs(f.values), [], 1);
    
    % Simplify (as suggested in Chebfun ticket #128)
    f = simplify(f);
    
    % Ensure f(-1) = 0:
    lval = get(f, 'lval');
    f.coeffs(end,:) = f.coeffs(end,:) - lval;
    f.values = bsxfun(@minus, f.values, lval);

end

function f = cumsumFiniteDim(f, m)
% CUMSUM over the finite dimension.

    for k = 1:m
        f.values = cumsum(f.values, 2);
        f.coeffs = cumsum(f.coeffs, 2);
        vscale = max(abs(f.values), [], 1);
        f.epslevel = sum(f.epslevel.*f.vscale, 2)/sum(vscale, 2); % TODO: Is this right?
        f.vscale = vscale;
    end
    out = f;

end


