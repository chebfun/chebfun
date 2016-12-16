function S = convertmat(n, K1, K2)
%CONVERTMAT  Conversion matrix used in the ultraspherical spectral method.
%   S = CONVERTMAT(N, K1, K2) computes the N-by-N matrix realization of the
%   conversion operator between two bases of ultrapherical polynomials.  The
%   matrix S maps N coefficients in a C^{(K1)} basis to N coefficients in a
%   C^{(K2 + 1)} basis, where, C^{(K)} denotes ultraspherical polynomial basis
%   with parameter K.  If K2 < K1, S is the N-by-N identity matrix.
%
%   This function is meant for internal use only and does not validate its
%   inputs.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Create the conversion matrix.
S = speye(n);
for s = K1:K2
    S = ultraS.spconvert(n, s) * S;
end

end
