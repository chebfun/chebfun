function values = coeffs2vals(coeffs)
%COEFFS2VALS   Convert Chebyshev coefficients to values at Chebyshev points
%of the 2nd kind.
%   V = COEFFS2VALS(C) returns the values of the polynomial V(i,1) = P(x_i) =
%   C(1,1)*T_{0}(x_i) + ... + C(N,1)*T_{N-1}(x_i), where the x_i are
%   2nd-kind Chebyshev nodes.
%
%   If the input C is an (N+1)xM matrix then V = COEFFS2VALS(C) returns the
%   (N+1)xM matrix of values V such that V(i,j) = P_j(x_i) = C(1,j)*T_{0}(x_i)
%   + C(2,j)*T_{1}(x_i) + ... + C(N,j)*T_{N-1}(x_i).
%
% See also VALS2COEFFS, CHEBPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Note]: This is equivalent to Discrete Cosine Transform of Type I.
%
% [Mathematical reference]: Sections 4.7 and 6.3 Mason & Handscomb, "Chebyshev
% Polynomials". Chapman & Hall/CRC (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *Note about symmetries* The code below takes steps to 
% ensure that the following symmetries are enforced:
% even Chebyshev COEFFS exactly zero ==> VALUES are exactly odd
% odd Chebychev COEFFS exactly zero ==> VALUES are exactly even
% These corrections are required because the MATLAB FFT does not
% guarantee that these symmetries are enforced.

% Get the length of the input:
n = size(coeffs, 1);

% Trivial case (constant or empty):
if ( n <= 1 )
    values = coeffs; 
    return
end

% check for symmetry
isEven = max(abs(coeffs(2:2:end,:)),[],1) == 0;
isOdd = max(abs(coeffs(1:2:end,:)),[],1) == 0;

% Scale them by 1/2:
coeffs(2:n-1,:) = coeffs(2:n-1,:)/2;

% Mirror the coefficients (to fake a DCT using an FFT):
tmp = [ coeffs ; coeffs(n-1:-1:2,:) ];

if ( isreal(coeffs) )
    % Real-valued case:
    values = real(fft(tmp));
elseif ( isreal(1i*coeffs) )
    % Imaginary-valued case:
    values = 1i*real(fft(imag(tmp)));
else
    % General case:
    values = fft(tmp);
end

% Flip and truncate:
values = values(n:-1:1,:);

% enforce symmetry
values(:,isEven) = (values(:,isEven)+flipud(values(:,isEven)))/2;
values(:,isOdd) = (values(:,isOdd)-flipud(values(:,isOdd)))/2;

end
