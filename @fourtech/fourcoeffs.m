function out = fourcoeffs(f, N)
%FOURCOEFFS   Fourier coefficients of a FOURTECH.
%   C = FOURCOEFFS(F) returns a column vector with the Fourier
%   coefficients of F using complex-exponential form. Specifically for
%   N = length(F):
%   If N is odd
%       F(x) = C(1)*z^(-(N-1)/2) + C(2)*z^(-(N-1)/2+1) + ... + C((N+1)/2) + ... 
%                + C(N-1)*z^((N-1)/2-1) + C(N)*z^((N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(-N/2) + C(2)*z^(-N/2+1) + ... + C(N/2+1) + ...
%                + C(N)*z^(N/2-1) + 
%   where z = exp(1i*pi*x).
%
%   A = FOURCOEFFS(F, N) truncates or pads the vector C so that N coefficients
%   of the FOURTECH F are returned.
%
%   If F is array-valued with M columns, then C is an NxM matrix.
%
% See also CHEBCOEFFS, POLY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    N = length(f);
end

if ( isempty(N) || (N <= 0) )
    out = [];
    return
end

% Get the coefficients:
c = f.coeffs;
[numCoeffs, numCols] = size(c);
fIsEven = ( mod(numCoeffs, 2) == 0 );
NisEven = ( mod(N, 2) == 0 ) ;

% Pad the coefficient vector:
if ( numCoeffs < N )
    c = [zeros(ceil((N-numCoeffs)/2), numCols) ; 
         c ; 
         zeros(ceil((N-numCoeffs)/2), numCols)];
end
numCoeffs = size(c, 1);

% Determine which index corresponds to the constant term:
if ( fIsEven )
    constIndex = numCoeffs/2;
else
    constIndex = (numCoeffs+1)/2;
end

% Use symmetry:
if ( NisEven )
    id = (constIndex-(N/2-1)) : (constIndex+(N/2));
else
    id = (constIndex-((N-1)/2)) : (constIndex+((N-1)/2));
end
c = c(id,:);

% Extract out the entries and flip the result to match the ordering from 
% negative powers to positive powers.
out = flipud(c);

end
