function values = chebpolyval(coeffs)
%CHEBPOLYVAL   Convert Chebyshev coefficients to values at Chebyshev points
%of the 1st kind.
%
%   V = CHEBPOLYVAL(C) returns the values of the polynomial V(i,1) = P(x_i) =
%   C(1,1)*T_{N-1}(x_i) + C(2,1)*T_{N-2}(x_i) + ... + C(N,1), where the x_i are
%   1st-kind Chebyshev nodes.
%
%   If the input C is an (N+1)*M matrix then V = CHEBPOLYVAL(C) returns the
%   (N+1)xM matrix of values V such that V(i,j) = P_j(x_i) = C(1,j)*T_{N-1}(x_i)
%   + C(2,j)*T_{N-2}(x_i) + ... + C(N,j)
%
%   See also CHEBPOLY, CHEBPTS.

%   Copyright 2013 by The University of Oxford and The Chebfun Developers. 
%   See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%   [Mathematical reference]

% Get the length of the input
n = size(coeffs, 1);

% Trivial case (constant)
if ( n == 1 )
    values = coeffs; 
    return
end

% Flip the input data up to down. The first row of the resulting data 
% should correspond to rightmost Chebyshev point of the 1st kind in [-1 1]. 

if ( isreal(coeffs) )
    values = realcoeffs(coeffs);
elseif ( isreal(1i*coeffs) )
    values = 1i*realcoeffs(imag(coeffs));
else
    values = realcoeffs(real(coeffs)) + 1i*realcoeffs(imag(coeffs));
end

end

% Real case - Chebyshev points of the 1st kind:
function v = realcoeffs(c)
n = size(c, 1);
m = size(c, 2);

coeffs = flipud(c); 
w = exp(1i*(0:n-1)*pi/(2*n)).';
if ( m > 1 )    % [TODO:] Delete this if block
    w = repmat(w, 1, m);
end
coeffs = w.*coeffs;
vv = n*real(ifft(coeffs));
if ( rem(n, 2) == 0 ) % Even case
    v(1:2:n-1,:) = vv(1:n/2,:);
    v(n:-2:2,:) = vv(n/2+1:n,:);
else                  % Odd case
    v(1:2:n,:) = vv(1:(n+1)/2,:);
    v(n-1:-2:2,:) = vv((n+1)/2+1:n,:);
end
v = flipud(v);

end
