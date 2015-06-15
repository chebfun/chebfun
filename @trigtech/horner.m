function y = horner(x, c, isReal)
%HORNER   Horner's scheme for evaluating a trigonometric expansion.
%   If C is a N x 1 column vector containing the coefficients for a trigonometric
%   series in complex exponential form, Y = HORNER(X, C) evaluates the 
%   trigonometric expansion
%
%     Y(x) = P_N(x) = C(1)*z^(-(N-1)/2) + C(2)*z^(-(N-1)/2-1) + ... + C(N)*z^((N-1)/2)
%   
%   if N is odd and 
%
%     Y(x) = P_N(x) = C(1)*z^(-N/2) + C(2)*z^(-N/2+1) + ... + C(N)*z^(N/2-1) 
%  
%   if N is even using Horner's method.  Here z = exp(1i*pi*x).
%
%   If C is an N x M matrix, then HORNER interprets each of the columns
%   of C as coefficients of a trigonometric series and evaluates the M 
%   Trigonmetric expansions at x.  The result is columns of a matrix Y =
%   [Y_1 ... Y_M].
%
%   Y = HORNER(X, C, ISREAL) uses a Horner scheme in real arithmetic if
%   ISREAL is true for all columns of C.
%
%   In both cases, X must be a column vector.
%
% See also TRIGTECH.FEVAL, TRIGTECH.BARY.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: horners is not typically called directly, but by FEVAL().
%
% Developer note: The algorithm is implemented both for scalar and for vector
% inputs and for both for cases where the expansions are known to be for a
% real-valued and complex-valued function. The vector and complex
% implementation could also be used for the scalar case, but the additional
% overheads make it a factor of 2-4 slower. Since the code is short, we
% live with the minor duplication.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    isReal = false;
end

% X should be a column vector.
if ( size(x, 2) > 1 )
    warning('CHEBFUN:TRIGTECH:horner:xDim', ...
        'Evaluation points should be a column vector.');
    x = x(:);
end

scalarValued = size(c,2) == 1;

if scalarValued && all(isReal)
    y = horner_scl_real(x, c);
elseif ~scalarValued && all(isReal)
    y = horner_vec_real(x, c);
elseif scalarValued && ~all(isReal)
    y = horner_scl_cmplx(x, c);
else
    y = horner_vec_cmplx(x, c);
    y(:,isReal) = real(y(:,isReal));
end

end

% Scalar version of Horner's scheme for evaluating one polynomial at the
% same locations.  Turns out to measurably faster than calling
% horner_vec_cmplx(x, c).
function q = horner_scl_cmplx(x, c)
N = size(c, 1);
z = exp(1i*pi*x); % Use complex exponential
q = c(N);  % Same as polyval for coeffs flipped from 0th deg to Nth deg
for j = N-1:-1:2
    q = z.*q + c(j);
end
if mod(N, 2) == 1
    q = exp(-1i*pi*(N-1)/2*x).*(z.*q + c(1));
else
    q = exp(-1i*pi*(N/2-1)*x).*q + cos(N/2*pi*x)*c(1);
end

end

% Vectorized version of Horner's scheme for evaluating multiple
% polynomials of the same degree at the same locations.
function q = horner_vec_cmplx(x, c)
nValsX = size(x, 1);
N = size(c, 1);
z = exp(1i*pi*x); % Use complex exponential
numCols = size(c,2);
if nValsX > 1
    z = repmat(z,[1 numCols]);
    e = ones(nValsX,1);
    q = repmat(c(N,:),[nValsX 1]);
    for j = N-1:-1:2
        q = z.*q + e*c(j,:);
    end
    if mod(N, 2) == 1
        q = bsxfun(@times, exp(-1i*pi*(N-1)/2*x), z.*q + e*c(1,:));
    else
        q = bsxfun(@times, exp(-1i*pi*(N/2-1)*x), q) + cos(N/2*pi*x)*c(1,:);
    end
else
    q = c(N,:);
    for j = N-1:-1:2
        q = z*q + c(j,:);
    end
    if mod(N, 2) == 1
        q = exp(-1i*pi*(N-1)/2*x)*(z.*q + c(1,:));
    else
        q = exp(-1i*pi*(N/2-1)*x)*q + cos(N/2*pi*z)*c(1,:);
    end
end

end

% Uses real arithmetic by summing the cosine/sine series version of the
% complex exponential series.  This is related to Clenshaw summation, but
% is not always called this in the trigonometric series context.
function q = horner_scl_real(x, c)

N = size(c, 1);

% Get all negative indexed coefficients coefficients so that the constant
% is the first term.
n = ceil((N+1)/2);
c = c(n:-1:1,:);
a = real(c);
b = imag(c);

% Just return the constant term.
if N == 1
    q = a;
    return
end
    
u = cos(pi*x);
v = sin(pi*x);
co = a(n);
si = b(n);
for j = n-1:-1:2
    temp = u.*co + v.*si + a(j);
    si = u.*si - v.*co + b(j);
    co = temp;
end

q = 2*u.*co + 2*v.*si + a(1);

end

% Uses real arithmetic by summing the cosine/sine series version of the
% complex exponential series.  This is related to Clenshaw summation, but
% is not always called this in the trigonometric series context.
function q = horner_vec_real(x, c)

nValsX = size(x, 1);
N = size(c, 1);
numCols = size(c,2);

% Get all negative indexed coefficients coefficients so that the constant
% is the first term.
n = ceil((N+1)/2);
c = c(n:-1:1,:);
a = real(c);
b = imag(c);
e = ones(nValsX,1);

% Just return the constant term.
if N == 1
    q = e*a;
    return
end

u = repmat(cos(pi*x),[1 numCols]);
v = repmat(sin(pi*x),[1 numCols]);
co = e*a(n,:);
si = e*b(n,:);
for j = n-1:-1:2
    temp = u.*co + v.*si + e*a(j,:);
    si = u.*si - v.*co + e*b(j,:);
    co = temp;
end

q = 2*u.*co + 2*v.*si + e*a(1,:);

end