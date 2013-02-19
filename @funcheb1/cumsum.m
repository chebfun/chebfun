function f = cumsum(f)
%CUMSUM    Indefinite integral of a FUNCHEB1
% CUMSUM(F) is the indefinite integral of the FUNCHEB1 F.

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

coeffs = f.coeffs;                     % Obtain Cheb coeffs {c_r}
[n, m] = size(coeffs);
coeffs = [zeros(2,m) ; coeffs];        % Pad with zeros
newcoeffs = zeros(n-1, m);             % Initialize vector {C_r}
newcoeffs(1:n-1,:) = ...               % Compute C_(n+1) ... C_2
    (coeffs(3:end-1,:) - coeffs(1:end-3,:)) ./ repmat(2*(n:-1:2)', 1, m);
newcoeffs(n,:) = coeffs(end,:) - coeffs(end-2,:)/2;   % Compute C_1
v = ones(1,n); 
v(end-1:-2:1) = -1;
newcoeffs(n+1,:) = v*newcoeffs;                       % Compute C_0

% [TODO] Should we still do this? Perhaps SIMPLIFY?
% % Trim small coeffs, as suggested in #128
% tol = chebfunpref('eps')/norm(newcoeffs,inf);
% idx = find(abs(newcoeffs)>tol,1);
% newcoeffs(1:idx-1) = [];

% Recover values and attach to output.
f.values = funcheb1.chebpolyval(newcoeffs);
f.coeffs = newcoeffs;

% Update vscale.
f.vscale = max(f.vscale, norm(f.values,inf));

end
