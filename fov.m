function [f, lineSegs] = fov(A, pref)
%FOV   Field of values (numerical range) of matrix A.
%   F = FOV(A), where A is a square matrix, returns a CHEBFUN F with domain [0
%   2*pi]. The image F([0 pi]) will be a curve describing the boundary of the
%   field of values A, a convex region in the complex plane. If A is Hermitian,
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
%   be used in constructing F. Note that PREF.splitting will always be set to
%   TRUE by FOV.
%
%   [F, LINESEGS] = FOV(A) also returns a cell array of structs defining the
%   line segments connecting up the extreme points: LINESEGS{k}.f is a CHEBFUN
%   with domain [-1, 1] connecting up the discontinuity in F at
%   LINESEGS{k}.theta.
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

if ( nargout == 1 )
    return
end

%% LINESEGS - implemented by Michael Overton
% Look for possible discontinuities at any break points. The boundary of the
% field of values may have corners, but it is continuous: hence the need
% sometimes for linear interpolation connecting up the sets of extreme points.

ends = f.domain;      % Break points are between 0 and 2*pi.
delta = 1e-14;        % Must be bigger than machine precision but not much 
                      % bigger: if it is too small spurious line segments will
                      % have tiny length, which hardly matters.        
if ( any(diff(ends) < delta) )
    delta = min(diff(ends))/3;
end
m = length(ends) - 1; % Number of pieces

% Get the values to the left and right of break points. The first break point is
% zero, which wraps around to 2*pi and needs special treatment. Exclude the
% right end point 2*pi.
lVal = feval(f, ends(1:end-1), 'left');  % Values to left of break points.
lVal(1) = feval(f, ends(end), 'left');   % Value to left of 2*pi.
rVal = feval(f, ends(1:end-1), 'right'); % Values to right of break points.

% Determine which points correspond to dicontinuities:
tol = 10*delta*max(abs([lVal, rVal]), [], 2);
discont = abs(lVal - rVal) > tol;

% Define additional CHEBFUNs for the line segments joining discontinuities.
% (First idea was to do this by inserting tiny intervals with steep slopes into
% the chebfun f, but this leads to loss of accuracy.)
indx = 0;
lineSegs = {};
for j = 1:m  
    if ( discont(j) )
        indx = indx + 1;
        % Line in the complex plane connecting the points lVal(j) and rVal(j):
        lineSegs{indx}.f = chebfun([lVal(j); rVal(j)]); % Default domain [-1 1].
        % The corresponding theta value in the Johnson algorithm:
        lineSegs{indx}.theta = ends(j); 
    end
end

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
