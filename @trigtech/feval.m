function y = feval(f, x)
%FEVAL   Evaluate a TRIGTECH.
%   Y = FEVAL(F, X) Evaluation of the TRIGTECH F at points X via
%   Horner's scheme.
%
%   If size(F, 2) > 1 then FEVAL returns values in the form [F_1(X), F_2(X),
%   ...], where size(F_k(X)) = size(X).
%
%   Example:
%     f = trigtech(@(x) exp(cos(pi*x)) );
%     x = linspace(-1, 1, 1000);
%     fx = feval(f, x);
%     plot(x,fx,'r-')
%

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    y = [];
    return 
end

% Reshape x to be a column vector for passing to polyval:
[~, m] = size(f);
fourierCoeff = f.coeffs(end:-1:1, :);
N = size(fourierCoeff, 1);
sizex = size(x);
ndimsx = ndims(x);
x = x(:);

if m > 1
    % If N is odd, no modification is needed.
    if ( mod(N, 2) == 1 )
       % Compute using vectorized polyval (Horner's scheme):
       y = bsxfun(@times, exp(-1i*pi*(N-1)/2*x), polyvalv(fourierCoeff, exp(1i*pi*x)));
    else % The degree N/2 term needs to be handled separately.   
       % Compute using vectorized polyval (Horner's scheme):
       y = cos(N/2*pi*x)*fourierCoeff(N,:) + ...
           bsxfun(@times, exp(-1i*pi*(N/2-1)*x), polyvalv(fourierCoeff(1:N-1,:), exp(1i*pi*x)));
    end
else
    % If N is odd, no modification is needed.
    if ( mod(N, 2) == 1 )
       % Compute using polyval (Horner's scheme):
       y = exp(-1i*pi*(N-1)/2*x).*polyval(fourierCoeff, exp(1i*pi*x));
    else % The degree (N/2 term) needs to be handled separately.
       % Compute using polyval (Horner's scheme):
       y = cos(N/2*pi*x)*fourierCoeff(N) + ...
           exp(-1i*pi*(N/2-1)*x).*polyval(fourierCoeff(1:N-1), exp(1i*pi*x));
    end
end    

% Return a real result wherever f is real.
y(:,f.isReal) = real(y(:,f.isReal));

% Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex(1), m*numel(x)/sizex(1));
end

% Vectorized version of Horner's scheme for evaluating multiple
% polynomials of the same degree at the same locations.
function q = polyvalv(c, x)
nValsX = size(x, 1);
degreePoly = size(c, 1);
q = ones(nValsX, 1)*c(1,:);
for j = 2:degreePoly
    q = bsxfun(@plus, bsxfun(@times, x, q), c(j,:));
end
end

end
