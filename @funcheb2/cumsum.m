function f = cumsum(f, pref)
%CUMSUM    Indefinite integral of a FUNCHEB2
%   CUMSUM(F) is the indefinite integral of the FUNCHEB2 F. The unknown constant
%   of integration is chosen so that F(-1) = 0.
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

% Trivial case of an empty funcheb2:
if ( isempty(f) )
    return
end

% Initialise storage:
c = f.coeffs;                     % Obtain Cheb coeffs c = {c_r}
[n, m] = size(c);
c = [zeros(2,m) ; c];             % Pad with zeros
C = zeros(n-1, m);                % Initialize vector C = {C_r}

% Compute C_(n+1) ... C_2:
C(1:n-1,:) = (c(3:end-1,:) - c(1:end-3,:)) ./ repmat(2*(n:-1:2)', 1, m);
C(n,:) = c(end,:) - c(end-2,:)/2; % Compute C_1
v = ones(1,n); 
v(end-1:-2:1) = -1;
C(n+1,:) = v*C;                   % Compute C_0

% Recover values and attach to output:
f.values = funcheb2.chebpolyval(C);
f.coeffs = C;

% Update vscale:
f.vscale = max(f.vscale, norm(f.values, inf));

% Simplify! (as suggested in #128)
if ( nargin  == 1 )
    pref = funcheb2.pref;
end
f = simplify(f, pref);

% Ensure f(-1) = 0:
f.coeffs(1,:) = f.coeffs(1,:) - f.values(1,:);
f.values = bsxfun(@minus, f.values, f.values(1,:));

end
