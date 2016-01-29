function v = dlt(c)
%DLT  Discrete Legendre transform. 
%   V = DLT(C) returns a column vector V such that
%       V(k) = C(1)*P_0(x(k)) + C(2)*P_1(x(k)) + ... C(N)*P_{N-1}(x(k)), 
%   where P_j(x) is the degree j Legendre polynomial and x is the vector of
%   Gauss-Legendre nodes (as returned by x = legpts(size(C,1))).
%
% See also CHEBFUN.IDLT, CHEBFUN.DCT.

% Nick Hale & Alex Townsend, Feb 2015.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: See N Hale & A Townsend - "A fast FFT-based DLT" for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(c, 1);

if ( N == 0 )
    
    % Trivial empty case.
    v = c;
    
elseif ( N < 5000 ) % <-- Determined expoerimentally.
    
    % Compute using direct (recurrence relation-based) method:
    v = dlt_direct(c);
    
else
    
    % Stage 1: Convert from Legendre to Chebyshev.
    c_cheb = leg2cheb(c);
    % Stage 2: Compute non-uniform DCT.
    v = chebfun.ndct(c_cheb);
    
end

end

function v = dlt_direct(c)
%DLT_DIRECT  Evaluate Legendre-Vandermonde matrix times a vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: Uses recurrence relation but does not form the matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(c, 1);
% Trivial case:
if ( N == 1 ), v = 1 + 0*c; return, end

% Compute the Legendre nodes:
x = legpts(N);

% Support for when c is a matrix:
v = repmat(c(1,:), length(x), 1) + bsxfun(@times, c(2,:), x);

Pm1 = 1 + 0*x; P = x;   % P_0 and P_1.
for n = 1:(N-2) % Recurrence relation.
    Pp1 = (2-1/(n+1))*(P.*x) - (1-1/(n+1))*Pm1; 
    Pm1 = P;    P = Pp1;
    v = v + bsxfun(@times, c(n+2,:), P);
end

end
