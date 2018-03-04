function f = randnfunsphere(lambda, type)
%RANDNFUNSPHERE   Smooth random function on the unit sphere
%   F = RANDNFUNSPHERE(LAMBDA) returns a SPHEREFUN of maximum
%   wavelength of about 2pi/LAMBDA and standard normal distribution
%   N(0,1) at each point.  F is obtained from a combination of
%   spherical harmonics with random coefficients.
%
%   RANDNFUNSPHERE(LAMBDA, 'monochrome') is similar, but uses a
%   fixed-degree expansion so that all components have wave number
%   about equal to 2pi/LAMBDA. Note that it is possible to shorten
%   'monochrome' to just 'm'.
%
%   RANDNFUNSPHERE() uses the default value LAMBDA = 1.
%
% Examples:
%
%   f = randnfunsphere(0.2); std2(f), plot(f,'zebra')
%
%   f = randnfunsphere(0.2,'monochrome'); std2(f), plot(f,'zebra')
%
% See also RANDNFUN, RANDNFUN2, RANDNFUNDISK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Parse the inputs
if nargin == 0
    lambda = 1;
    type = 'white';
elseif nargin == 1
    % User either called randnspherefun(lambda) or randnspherefun('monochromatic')
    if isnumeric( lambda )
        type = 'white';
    elseif ( ischar( lambda ) && strncmpi(lambda,'m',1) )
        lambda = 1;
        type = 'monochromatic';
    else
        error('CHEBFUN:SPHEREFUN:randnspherefun:inputUnknown', ...
            'Input argument unknown, options are a positive number or the string ''monochromatic''.');
    end
elseif nargin == 2
    if ~isnumeric( lambda )
        error('CHEBFUN:SPHEREFUN:randnspherefun:inputUnknown', ...
            'First input argument must be a positive number.');
    end
    
    if ( ~ischar( type ) || ~strncmpi(type,'m',1) )
        error('CHEBFUN:SPHEREFUN:randnspherefun:inputUnknown', ...
            'Second input argument unknown.  The only option is ''monochromatic''.');
    else
        type = 'monochromatic';
    end
end

% Since the unit sphere has circumference L = 2*pi, the following
% choice matches the choice deg = floor(L/lambda) in randnfun.
deg = floor(2*pi/lambda);

% We do not use adaptive construction, but just sample the function on a
% fine enough grid to exactly resolve it then pass this to the constructor.

% Sampling grid to exactly recover the random spherical harmonic:
ll = trigpts(2*deg+2,[-pi pi]);
tt = linspace(0,pi,2*deg+2);

if strcmpi(type,'monochromatic');
    c = randn(2*deg+1, 1);
    c = sqrt(4*pi/nnz(c))*c;     % normalize so variance is 1
    f = spherefun( sphHarmSumFixedDeg(ll,tt,deg,c) );
else
    c = randn((deg+1)^2, 1);
    c = sqrt(4*pi/nnz(c))*c;         % normalize so variance is 1
    f = spherefun( sphHarmSum(ll,tt,deg,c) );
end

% Simplify the result
f = simplify(f);

end

function F = sphHarmSum(lam,th,deg,coeffs)
%SPHHARMSUM Linear combination of all spherical harmonics of a given degree 
%   
%   F = SPHHARMSUM(LAM,TH,DEG,COEFFS) computes the linear combination of
%   all spherical harmonics up to degree DEG over a tensor product grid
%   given by LAM x TH using the coefficients specified by COEFFS. LAM and
%   TH are assumed to be vectors containing slices from the tensor product
%   grid. COEFFS must contain (DEG+1)^2 entries, as this is the total
%   number of linearly independent spherical harmonics <= deg.  The
%   coefficients are assumed to be ordered corresponding to the degree and
%   order of the spherical harmonics.  For example, for a degree 3 sum, the
%   ordering should be 0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3.

% Make th a row vector to better work with Matlab's Legendre function. Also
% Legendre operates on cos(th).
costh = cos(th(:)).';
% Make lam row
lam = lam(:).';

% Handle the zero degree term separately since it's simple
c = coeffs(1);
coeffs = coeffs(2:end);
F = 1/sqrt(4*pi)*c;

for l = 1:deg
    % Normalization terms for the associated Legendre functions.
    m = (0:l).';
    a = (-1).^m./sqrt((1+double(m==0))*pi);

    % Compute the associated Legendre functions of cos(th) (Co-latitude)
    G = legendre(l, costh, 'norm');
    
    % Extract out the coefficients and then discard the coefficients from
    % the global coeffs vector.
    c = coeffs(1:(2*l+1));
    coeffs = coeffs((2*l+1)+1:end);

    % Multiply the random coefficients by the associated Legendre polynomials.
    % We will do one for the positive (including zero) order associated
    % Legendre functions and one for the negative.
    Gp = bsxfun(@times,a.*c(l+1:end),G);
    Gn = bsxfun(@times,a(2:end).*c(l:-1:1), G(2:end,:,:));
    
    % Permute rows and columns to do the tensor product computation
    Gp = permute(Gp,[2 3 1]);
    Gn = permute(Gn,[2 3 1]);
    
    % Multiply the associated Legendre polynomials by the correct Fourier modes
    % in the longitude variable and sum up the results.
    F = F + sum(bsxfun(@mtimes,Gp, permute(cos((0:l)'*lam),[3 2 1])),3) + ...
            sum(bsxfun(@mtimes,Gn, permute(sin((1:l)'*lam),[3 2 1])),3);
end

end

function F = sphHarmSumFixedDeg(lam,th,l,c)
%SPHHARMSUMFIXEDDEG Linear combination of all spherical harmonics of a fixed degree 
%   
%   F = SPHHARMSUMFIXEDDEG(LAM,TH,DEG,COEFFS) computes the linear
%   combination of all spherical harmonics of degree DEG over a tensor
%   product grid given by LAM x TH using the coefficients specified by
%   COEFFS. LAM and TH are assumed to be vectors containing slices from
%   the tensor product grid. COEFFS must contain 2*DEG+1 entries, as this
%   is the total number of linearly independent spherical harmonics of
%   degree DEG. The coefficients are assumed to be ordered corresponding to
%   order of the spherical harmonics.  For example, for a degree 3 sum, the
%   ordering should be -3,-2,-1,0,1,2,3.

% Make th a row vector to better work with Matlab's Legendre function. Also
% Legendre operates on cos(th).
costh = cos(th(:)).';
% Make lam row
lam = lam(:).';

% Normalization terms for the associated Legendre functions.
m = (0:l).';
a = (-1).^m./sqrt((1+double(m==0))*pi);

% Compute the associated Legendre functions of cos(th) (Co-latitude)
F = legendre(l, costh, 'norm');

% Multiply the random coefficients by the associated Legendre polynomials.
% We will do one for the positive (including zero) order associated
% Legendre functions and one for the negative.
Fp = bsxfun(@times,a.*c(l+1:end),F);
Fn = bsxfun(@times,a(2:end).*c(l:-1:1), F(2:end,:,:));

% Permute rows and columns to do the tensor product computation
Fp = permute(Fp,[2 3 1]);
Fn = permute(Fn,[2 3 1]);

% Multiply the associated Legendre polynomials by the correct Fourier modes
% in the longitude variable and sum up the results.
F = sum(bsxfun(@mtimes,Fp, permute(cos((0:l)'*lam),[3 2 1])),3) + ...
    sum(bsxfun(@mtimes,Fn, permute(sin((1:l)'*lam),[3 2 1])),3);

end

