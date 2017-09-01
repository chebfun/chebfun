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

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
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

% Just return the constant term.
if N == 1
    q = ones( size(x, 1), 1)*c(N);
    return
end

z = exp(1i*pi*x); % Use complex exponential
q = c(N);  % Same as polyval for coeffs flipped from 0th deg to Nth deg

for j = N-1:-1:2
    q = c(j) + z.*q;
end
if mod(N, 2) == 1
    q = exp(-1i*pi*(N-1)/2*x).*(c(1) + z.*q);
else
    q = exp(-1i*pi*(N/2-1)*x).*q + cos(N/2*pi*x)*c(1);
end

end

% Vectorized version of Horner's scheme for evaluating multiple
% polynomials of the same degree at the same locations.
function q = horner_vec_cmplx(x, c)
nValsX = size(x, 1);
N = size(c, 1);
% Just return the constant term.
if N == 1
    q = repmat(c(N,:),[nValsX 1]);
    return
end

z = exp(1i*pi*x); % Use complex exponential
numCols = size(c,2);
if nValsX > 1
    z = repmat(z,[1 numCols]);
    e = ones(nValsX,1);
    q = repmat(c(N,:),[nValsX 1]);
    for j = N-1:-1:2
        q = e*c(j,:) + z.*q;
    end
    if mod(N, 2) == 1
        q = bsxfun(@times, exp(-1i*pi*(N-1)/2*x), e*c(1,:) + z.*q);
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
        q = exp(-1i*pi*(N/2-1)*x)*q + cos(N/2*pi*x)*c(1,:);
    end
end

end

% Uses real arithmetic by summing the cosine/sine series version of the
% complex exponential series.
function q = horner_scl_real(x, c)

N = size(c, 1);

% Get all negative indexed coefficients coefficients so that the constant
% is the first term.
n = ceil((N+1)/2);
c = c(n:-1:1,:);
a = real(c);
b = imag(c);

% Adjust the last coefficient which corresponds to the pure cos(pi*n*x) 
% mode in the case that N is even.
if mod(N,2) == 0
    a(n) = a(n)/2;
    b(n) = 0;
end

% Just return the constant term.
if N == 1
    q = ones( size(x, 1), 1)*a;
    return
end
    
u = cos(pi*x);
v = sin(pi*x);
co = a(n);
si = b(n);
for j = n-1:-1:2
    temp = a(j) + u.*co + v.*si;
    si = b(j) + u.*si - v.*co;
    co = temp;
end

q = a(1) + 2*(u.*co + v.*si);

end

% Uses real arithmetic by summing the cosine/sine series version of the
% complex exponential series.  This is measurably faster when c is a
% matrix, i.e array valued polynomial.
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

% Adjust the last coefficient which corresponds to the pure cos(pi*n*x) 
% mode in the case that N is even.
if mod(N,2) == 0
    a(n,:) = a(n,:)/2;
    b(n,:) = 0;
end

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
    temp = e*a(j,:) + u.*co + v.*si;
    si = e*b(j,:) + u.*si - v.*co;
    co = temp;
end

q = e*a(1,:) + 2*(u.*co + v.*si);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The code below implements the Clenshaw algortithm for
% summing a Fourier series (also called the Goertzel-Watt algorithm).  This
% algorithm is about 1.5 times faster than the horner schemes implemented
% above.  However, it is unstable when x is close to zero or +/- 1 (and the
% periodic extensions of these values); see, for example, W. M. Gentleman,
% "An error analysis of Goertzel's (Watt's) method for computing Fourier
% coefficients" or A. C. R. Newbery, "Error analysis for Fourier series
% evaluation." Mathematics of Computation 27.123 (1973): 639-644.  The
% latter article describes how to fix this using a "phase-shift" for select
% values of x in the unstable region. However, the logic involved in doing
% this will probably outweigh any benefits of the speed-up of this
% algorithm over the Horner schemes implemented above.  I am leaving these
% algorithms here as someone in the future may figure out a slick way to
% implement stable versions of them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function q = clenshaw_scl_real(x, c)
% 
% N = size(c, 1);
% 
% % Get all negative indexed coefficients coefficients so that the constant
% % is the first term.
% n = ceil((N+1)/2);
% c = c(n:-1:1,:);
% a = real(c);
% b = imag(c);
% 
% % Adjust the last coefficient which corresponds to the pure cos(pi*n*x) 
% % mode in the case that N is even.
% if mod(N,2) == 0
%     a(n) = a(n)/2;
%     b(n) = 0;
% end
% 
% % Just return the constant term.
% if N == 1
%     q = a;
%     return
% end
% 
% % Clenshaw scheme for scalar-valued functions.
% bk1u = 0*x; 
% bk2u = bk1u;
% bk1v = bk1u; 
% bk2v = bk1u;
% 
% u = 2*cos(pi*x);
% n = size(b,1)-1;
% for k = (n+1):-2:3
%     bk2u = a(k) + u.*bk1u - bk2u;
%     bk1u = a(k-1) + u.*bk2u - bk1u;
%     bk2v = b(k) + u.*bk1v - bk2v;
%     bk1v = b(k-1) + u.*bk2v - bk1v;
% end
% if ( mod(n, 2) )
%     [bk1u, bk2u] = deal(a(2) + u.*bk1u - bk2u, bk1u);
%     bk1v = b(2) + u.*bk1v - bk2v;
% end
% q = a(1) + u.*bk1u - 2*bk2u + (2*sin(pi*x)).*bk1v;
% end
% 
% % Uses real arithmetic by summing the cosine/sine series version of the
% % complex exponential series.  This is related to Clenshaw summation, but
% % is not always called this in the trigonometric series context.
% function q = clenshaw_vec_real(x, c)
% 
% nValsX = size(x, 1);
% N = size(c, 1);
% numCols = size(c,2);
% 
% % Get all negative indexed coefficients coefficients so that the constant
% % is the first term.
% n = ceil((N+1)/2);
% c = c(n:-1:1,:);
% a = real(c);
% b = imag(c);
% e = ones(nValsX,1);
% 
% % Adjust the last coefficient which corresponds to the pure cos(pi*n*x) 
% % mode in the case that N is even.
% if mod(N,2) == 0
%     a(n,:) = a(n,:)/2;
%     b(n,:) = 0;
% end
% 
% % Just return the constant term.
% if N == 1
%     q = e*a;
%     return
% end
% 
% x = repmat(x,[1 numCols]);
% 
% % Clenshaw scheme for scalar-valued functions.
% bk1u = zeros(nValsX,numCols);
% bk2u = bk1u;
% bk1v = bk1u; 
% bk2v = bk1u;
% 
% u = 2*cos(pi*x);
% n = size(b,1)-1;
% for k = (n+1):-2:3
%     bk2u = e*a(k,:) + u.*bk1u - bk2u;
%     bk1u = e*a(k-1,:) + u.*bk2u - bk1u;
%     bk2v = e*b(k,:) + u.*bk1v - bk2v;
%     bk1v = e*b(k-1,:) + u.*bk2v - bk1v;
% end
% if ( mod(n, 2) )
%     [bk1u, bk2u] = deal(e*a(2,:) + u.*bk1u - bk2u, bk1u);
%     bk1v = e*b(2,:) + u.*bk1v - bk2v;
% end
% q = e*a(1,:) + u.*bk1u - 2*bk2u + (2*sin(pi*x)).*bk1v;
% end
