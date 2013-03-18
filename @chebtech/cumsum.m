function f = cumsum(f, m, pref)
%CUMSUM    Indefinite integral of a CHEBTECH.
%   CUMSUM(F) is the indefinite integral of the CHEBTECH F with the constant of
%   integration is chosen so that F(-1) = 0. 
%
%   CUMSUM(F, M) will compute the Mth definite integral with the constant of
%   integration chosen so that each intermediary integral evaluates to 0 at -1
%   so that CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, PREF) or CUMSUM(F, M,  PREF) uses options from the preference
%   structure PREF when building the output CHEBTECH.
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the fun G of length n is represented as
%       SUM_{r=0}^{n-1} c_r T_r(x)
% its integral is represented with a fun of length n+1 given by
%       SUM_{r=0}^{n} C_r T_r (x)
% where C_0 is determined from the constant of integration as
%       C_0 = SUM_{r=1}^{n} (-1)^(r+1) C_r;
% C_1 = c_0 - c_2/2, and for r > 0,
%       C_r = (c_{r-1} - c_{r+1})/(2r),
% with c_{n+1} = c_{n+2} = 0.
% (See "Chebyshev Polynomials" by Mason and Handscomb, CRC 2002, pg 32-33.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trivial case of an empty CHEBTECH:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    m = 1;
    pref = chebtech.pref;
elseif ( nargin < 3 )
    if ( isstruct(m) )
        pref = m;
        m = 1;
    else
        pref = chebtech.pref();
    end
end

% Initialise storage:
c = f.coeffs;                         % Obtain Cheb coeffs c = {c_r}

% Loop for higher-order integrals:
for k = 1:m

    [n, m] = size(c);
    c = [zeros(2,m) ; c];             % Pad with zeros
    C = zeros(n-1, m);                % Initialize vector C = {C_r}

    % Compute C_(n+1) ... C_2:
    C(1:n-1,:) = (c(3:end-1,:) - c(1:end-3,:)) ./ repmat(2*(n:-1:2)', 1, m);
    C(n,:) = c(end,:) - c(end-2,:)/2; % Compute C_1
    v = ones(1,n);
    v(end-1:-2:1) = -1;
    C(n+1,:) = v*C;                   % Compute C_0 (satisfies f(-1) = 0)
    
    % Copy coeffs back into c for loop:
    c = C;
end

% Recover values and attach to output:
f.values = f.chebpolyval(C);
f.coeffs = C;

% Update vscale:
f.vscale = max(f.vscale, max(abs(f.values), [], 1));

% Simplify! (as suggested in #128)
if ( nargin  == 1 )
    pref = chebtech.pref();
end
f = simplify(f, pref);

% Ensure f(-1) = 0:
lval = get(f, 'lval');
f.coeffs(end,:) = f.coeffs(end,:) - lval;
f.values = bsxfun(@minus, f.values, lval);

end
