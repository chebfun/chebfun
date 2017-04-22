function f = randnfunsphere(dt)
%RANDNFUN   Random smooth function on the unit sphere
%   F = RANDNFUNSPHERE(DT) returns a smooth SPHEREFUN of maximum
%   frequency about 2pi/DT and standard normal distribution N(0,1)
%   at each point.  F is obtained from a combination of spherical
%   harmonics with random coefficients.
%
%   F = RANDNFUNSPHERE() uses the default value DT = 1.
%
% Examples:
%
%   f = randnfunsphere(0.2); std2(f), plot(f)
%   colormap([0 0 0; 1 1 1]); caxis(norm(caxis,inf)*[-1 1])
%
% See also RANDNFUN, RANDNFUN2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 0
    dt = 1;
end

deg = round(pi/dt);
c = randn(2*deg+1, 1);
c = 2*c/sqrt(nnz(c));     % normalize so variance is 1
f = spherefun(@(lam, th) sphHarmFixedDegRand(lam, th, deg, c));

end

function F = sphHarmFixedDegRand(lam,th,l,c)
%SPHHARMFIXEDDEGRAND Random combination of all spherical harmonics of fixed degree
%   F = SPHHARMFIXEDDEGRAND(LAM,TH,DEG,C) is a random combination of all
%   spherical harmonics of a given degree DEG.  The random coefficients
%   for the combination are given in C, which must be of dimension 2*DEG+1.

% Determine whether the input is on a tensor product grid.  If it is then
% we can speed things up because the associated Legendre functions can be
% computed more quickly.

[m, n] = size(th);
tensorGrid = 0;
if m > 1 && n > 1
    th = th(:,1);
    tensorGrid = 1;
else
    % Flatten theta so it works with Matlab's Legendre function
    th = th(:).';
end

% Flatten lambda so it works with Matlab's Legendre function
lam = lam(:).';

% Normalization terms for the associated Legendre functions.

mm = ones(l+1, 1)*(1:2*l);
pp = (0:l)'*ones(1, 2*l);
mask = (pp > abs(mm-(2*l+1)/2));
kk = mm.*mask + ~mask;
aa = exp(-sum(log(kk), 2));
a = 2*sqrt((2*l + 1)/4/pi*aa);
a(1) = a(1)/2;                      % correction for the zero mode.

% Compute the associated Legendre functions of cos(th) (Co-latitude)
F = legendre(l, cos(th));

% Multiply the random coefficients by the associated Legendre polynomials.
% We will do one for the positive (including zero) order associated
% Legendre functions and one for the negative.
Fp = bsxfun(@times,a.*c(l+1:end),F);
Fn = bsxfun(@times,a(2:end).*c(l:-1:1), F(2:end,:,:));

% If this is a tensor grid then reproduce the associated Legendre functions
% to match the tensor grid structure.
if ( tensorGrid )
    Fp = repmat(Fp, [1 n]);
    Fn = repmat(Fn, [1 n]);
end

% Multiply the associated Legendre polynomials by the correct Fourier modes
% in the longitude variable and sum up the results.
F = sum(Fp.*cos((0:l)'*lam)) + sum(Fn.*sin((1:l)'*lam));
% Reshape so it is the same size as the th and lam that were passed in.
F = reshape(F, [m n]);

end
