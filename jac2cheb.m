function c_cheb = jac2cheb(c_jac, alph, bet)
%JAC2CHEB convert Legendre coefficients to Chebyshev coefficients. 
%   C_CHEB = JAC2CHEB(C_JAC, A, B) converts the vector C_JAC of Jacobi
%   coefficients to a vector C_CHEB of Chebyshev coefficients such that
%       C_CHEB(1)*T_0(x) + ... + C_CHEB(N)*T{N-1}(x) = ...
%           C_JAC(1)*P_0^{(A,B)}(x) + ... + C_JAC(N)*P_{N-1}^{(A,B)},
%   where P_k^{(A,B)} is the degree k Jacobi polynomial corresponding to the
%   weight function w(x) = (1-X)^A * (1+X)^B.
%
%   If C_LEG is a matrix then the LEG2CHEB operation is applied to each column.
%
% See also CHEB2JAC. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  This simply uses the recurrence relation to for the Jacobi-Vandermode matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( alph == 0 && bet == 0 )
    % Special case:
    c_cheb = leg2cheb(c_jac); 
    return
end

c_cheb = jac2cheb_direct(c_jac, alph, bet); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECT METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c_cheb = jac2cheb_direct(c_jac, a, b)
%JAC2CHEB_DIRECT   Convert Leg to Cheb coeffs using the 3-term recurrence.
N = size(c_jac,1) - 1;                      % Degree of polynomial.
% Don't let N be too big:
if ( N > 2^11 )
    error('CHEBFUN:jac2cheb:arraySize', ...
        'Maximum transform size (2048) is exceeded.');
end
if ( N <= 0 ), c_cheb = c_jac; return, end  % Trivial case.
x = cos(pi*(0:N)'/N);                       % Chebyshev grid (reversed order).
% Make the Jacobi-Chebyshev Vandemonde matrix:
apb = a + b; aa  = a * a; bb  = b * b;
P = zeros(N+1); P(:,1) = 1;    
P(:,2) = 0.5*(2*(a + 1) + (apb + 2)*(x - 1));   
for k = 2:N
    k2 = 2*k;
    k2apb = k2 + apb;
    q1 =  k2*(k + apb)*(k2apb - 2);
    q2 = (k2apb - 1)*(aa - bb);
    q3 = (k2apb - 2)*(k2apb - 1)*k2apb;
    q4 =  2*(k + a - 1)*(k + b - 1)*k2apb;
    P(:,k+1) = ((q2 + q3*x).*P(:,k) - q4*P(:,k-1)) / q1;
end
v_cheb = P*c_jac;                           % Values on Chebyshev grid.
c_cheb = idct1(v_cheb);                     % Chebyshev coefficients.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DCT METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = idct1(v)
%IDCT1   Convert values on a Cheb grid to Cheb coefficients (inverse DCT1).
% IDCT1(V) returns T_N(X_N)\V, where X_N = cos(pi*(0:N))/N and T_N(X) = [T_0,
% T_1, ..., T_N](X) where T_k is the kth 1st-kind Chebyshev polynomial.
N = size(v, 1);                             % Number of terms.
c = fft([v ; v(N-1:-1:2,:)])/(2*N-2);       % Laurent fold and call FFT.
c = c(1:N,:);                               % Extract the first N terms.
if (N > 2), c(2:N-1,:) = 2*c(2:N-1,:); end  % Scale interior coefficients.
if ( isreal(v) ), c = real(c); end          % Ensure a real output.
end
