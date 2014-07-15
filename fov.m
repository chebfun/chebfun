function f = fov(A, pref)
%FOV   Field of values (numerical range) of matrix A.
%   F = FOV(A), where A is a square matrix, returns a CHEBFUN F with domain [0
%   2*pi].  The image F([0 pi]) will be a curve describing the boundary of the
%   field of values A, a convex region in the complex plane.  If A is Hermitian,
%   the field of values is a real interval, and if A is normal, it is the convex
%   hull of the eigenvalues of A.
%
%   The numerical abscissa of A is equal to max(real(F)), though this is much
%   better computed as max(real(eig(A + A')))/2.
%
%   The algorithm use is that of C. R. Johnson, Numerical determination of the
%   field of values of a general complex matrix, SIAM J. Numer. Anal. 15 (1978),
%   595-602.
%
%   F = FOV(A, PREF) allows the preferences in the CHEBFUNPREF structure PREF to
%   be used in constructing F. Note that PREF.splitting will
%   always be set to TRUE by FOV.
%
% Example:
%   A = randn(5);
%   F = fov(A);
%   hold off, fill(real(F), imag(F), [1 .5 .5]), axis equal
%   e = eig(A);
%   hold on, plot(real(e), imag(e), '.k', 'markersize', 16)

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    % Obtain preferences:
    pref = chebfunpref();
end

% Set splitting to 'on' and the domain to [0, 2*pi].
pref.splitting = true;
pref.domain = [0, 2*pi];

% Construct a CHEBFUN of the FOV curve:
f = chebfun(@(theta) fovCurve(theta, A), pref);

end

function z = fovCurve(theta, A)
z = NaN(size(theta));
for j = 1:length(theta)
    r = exp(1i*theta(j));
    B = r*A;
    H = (B + B')/2;
    [X, D] = eig(H);
    [ignored, k] = max(diag(D));
    v = X(:,k);
    z(j) = v'*A*v/(v'*v);
end
end
