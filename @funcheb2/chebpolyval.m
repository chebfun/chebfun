function values = chebpolyval(coeffs)
%CHEBPOLYVAL   Convert Chebyshev coefficients to values at Chebyshev points.
%   V = CHEBPOLYVAL(C) returns the values of the polynomial V(i,1) = P(x_i) =
%   C(1,1)*T_{N-1}(x_i) + C(2,1)*T_{N-2}(x_i) + ... + C(N,1), where the x_i are
%   2nd-kind Chebyshev nodes.
%
%   If the input C is an (N+1)*M matrix, then V = chebpolyval(C) returns the
%   (N+1)xM matrix of values V such that V(i,j) = P_j(x_i) = C(1,j)*T_{N-1}(x_i)
%   + C(2,j)*T_{N-2}(x_i) + ... + C(N,j)
%
%   See also chebpoly, chebpts.
%
%   [TODO]: [Mathmatical reference]

%   Copyright 2013 by The University of Oxford and The Chebfun Developers. 
%   See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Get the length of the input
n = size(coeffs, 1);

% Trivial case (constant)
if ( n == 1 )
    values = coeffs; 
    return
end

% Compute an index for and scale the interior coefficients
ii = 2:n-1;
coeffs(ii,:) = 0.5*coeffs(ii,:);

% Mirror the coefficients (to fake a DCT using an FFT)
tmp = [coeffs(end:-1:1,:) ; coeffs(ii,:)];

if ( isreal(coeffs) )
    % Real-valued case
    values = real(ifft(tmp));
elseif ( isreal(1i*coeffs) )
    % Imaginary-valued case
    values = 1i*real(ifft(imag(tmp)));
else
    % General case
    values = ifft(tmp);
end

% Scaling
values = (n-1)*[ 2*values(1,:) ; values(ii,:) + values(2*n-ii,:) ; 2*values(n,:) ];

% Flip the order
values = values(end:-1:1,:);

end

