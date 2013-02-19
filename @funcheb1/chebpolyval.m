function values = chebpolyval(coeffs)
%CHEBPOLYVAL   Convert Chebyshev coefficients to values at Chebyshev points
%   of 1st kind.
%
%   V = CHEBPOLYVAL(C) returns the values of the polynomial 
%          P_j(x) = C(1,j)*T_{N-1}(x) + C(2,j)*T_{N-2}(x) + ... + C(N,j) 
%   at 1st-kind Chebyshev nodes so that V(i,j) = P_j(x_i).
%
%   See also chebpoly, chebpts.
%
%   [Mathmatical reference]

%   Copyright 2013 by The University of Oxford and The Chebfun Developers. 
%   See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Get the length of the input
n = size(coeffs, 1);

% Trivial case (constant)
if ( n == 1 )
    values = coeffs; 
    return
end

% Flip the input data up to down. The first row of the resulting data 
% should correspond to rightmost Chebyshev point of 1st kind in [-1 1]. 

if isreal(coeffs)
    values = realcoeffs(coeffs);
elseif isreal(1i*coeffs)
    values = 1i*realcoeffs(imag(coeffs));
else
    values = realcoeffs(real(coeffs))+1i*realcoeffs(imag(coeffs));
end

end

function v = realcoeffs(c) % Real case - Chebyshev points of the 1st kind
n = size(c,1);
m = size(c,2);
%%%%%%%%% without FFT %%%%%%%%%%%%%%%%%%%
% c = flipud(c);
% k = 0:n-1;
% K = pi*k/(2*n);
% v = zeros(1,n);
% 
% for j = 0:n-1
%        v(j+1) = cos(K*(2*j+1))*c;    
% end
% v = v.';
% v = v(n:-1:1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %b = repmat(exp(-1i*(n-1:-1:0)*pi/(2*n)).',1,m);
% N = n-1;
% b = exp(-1i*pi*(0:N)'/(N+1));
% c = flipud(c);
% C = c.*b;
% C = [C;-C];
% 
% % Mirror the data
% %c = [c(n:-1:1,:); zeros(1,m); c(1:n-1,:)];
% % c(1,:) = c(1,:)/2; c(n,:) = c(n,:)/2;
% v = real(fft(C));
% v = v(n:-1:1,:)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coeffs = flipud(c); 
w = exp(1i*(0:n-1)*pi/(2*n)).';
if m > 1
    w = repmat(w, 1, m);
end
coeffs = w.*coeffs;
vv = n*real(ifft(coeffs));
if ( rem(n,2) == 0 ) % Even case
    v(1:2:n-1,:) = vv(1:n/2,:);
    v(n:-2:2,:) = vv(n/2+1:n,:);
else             % Odd case
    v(1:2:n,:) = vv(1:(n+1)/2,:);
    v(n-1:-2:2,:) = vv((n+1)/2+1:n,:);
end
v = flipud(v);

end
