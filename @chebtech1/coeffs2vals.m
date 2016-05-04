function values = coeffs2vals(coeffs)
%COEFFS2VALS   Convert Chebyshev coefficients to values at Chebyshev points
%of the 1st kind.
%   V = COEFFS2VALS(C) returns the values of the polynomial V(i,1) = P(x_i) =
%   C(1,1)*T_{N-1}(x_i) + C(2,1)*T_{N-2}(x_i) + ... + C(N,1), where the x_i are
%   1st-kind Chebyshev nodes.
%
%   If the input C is an (N+1)xM matrix then V = COEFFS2VALS(C) returns the
%   (N+1)xM matrix of values V such that V(i,j) = P_j(x_i) = C(1,j)*T_{N-1}(x_i)
%   + C(2,j)*T_{N-2}(x_i) + ... + C(N,j)
%
% See also VALS2COEFFS, CHEBPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Note]: This is equivalent to Discrete Cosine Transform of Type III.
%
% [Mathematical reference]: Section 4.7 Mason & Handscomb, "Chebyshev
% Polynomials". Chapman & Hall/CRC (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *Note about symmetries* The code below takes steps to 
% ensure that the following symmetries are enforced:
% even Chebyshev COEFFS exactly zero ==> VALUES are exactly odd
% odd Chebychev COEFFS exactly zero ==> VALUES are exactly even
% These corrections are required because the MATLAB FFT does not
% guarantee that these symmetries are enforced.

% Get the length of the input:
[n, m] = size(coeffs);

% Trivial case (constant):
if ( n <= 1 )
    values = coeffs;
    return
end

% check for symmetry
isEven = max(abs(coeffs(2:2:end,:)),[],1) == 0;
isOdd = max(abs(coeffs(1:2:end,:)),[],1) == 0;

% Computing the weight vector often accounts for at least half the cost of this
% transformation. Given that (a) the weight vector depends only on the length of
% the coefficients and not the coefficients themselves and (b) that we often
% perform repeated transforms of the same length, we store w persistently.
persistent w
if ( size(w, 1) ~= 2*n )
    % Pre-compute the weight vector:
    w = (exp(-1i*(0:2*n-1)*pi/(2*n))/2).';
    w([1, n+1]) = [2*w(1); 0];
    w(n+2:end) = -w(n+2:end);
end

% Mirror the values for FFT:
c_mirror = [coeffs ; ones(1, m) ; coeffs(end:-1:2,:)];

% Apply the weight vector:
c_weight = bsxfun(@times, c_mirror, w);
values = fft(c_weight);

% Truncate and flip the order:
values = values(n:-1:1,:);

% Post-process:
if ( isreal(coeffs) )           
    % Real-valued case:
    values = real(values);
elseif ( isreal(1i*coeffs) )    
    % Imaginary-valued case:
    values = 1i*imag(values);
end

% enforce symmetry
values(:,isEven) = (values(:,isEven)+flipud(values(:,isEven)))/2;
values(:,isOdd) = (values(:,isOdd)-flipud(values(:,isOdd)))/2;

end
