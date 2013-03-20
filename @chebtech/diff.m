function f = diff(f, k, dim)
%DIFF   Derivative of a CHEBTECH.
%   DIFF(F) is the derivative of F and DIFF(F, K) is the Kth derivative.
%
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference along
%   dimension DIM. For DIM = 1, this is the same as above. For DIM = 2, this
%   is a finite difference along the columns of a vectorised CHEBTECH.
%   If F has L columns, an empty CHEBTECH will be returned for K >= L.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the CHEBTECH G of length n is represented as
%       \sum_{r=0}^{n-1} C_r T_r(x)
% its derivative is represented with a CHEBTECH of length n-1 given by
%       \sum_{r=0}^{n-2} c_r T_r (x)
% where c_0 is determined by
%       c_0 = c_2/2 + C_1;
% and for r > 0,
%       c_r = c_{r+2} + 2*(r+1)*C_{r+1},
% with c_n = c_{n+1} = 0.
% (See "Chebyshev Polynomials" by Mason and Handscomb, CRC 2002, p. 34.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check the inputs:

% Trivial case of an empty CHEBTECH:
if ( isempty(f) )
    return
end

if ( nargin < 2 || isempty(k) )
    % Assume 1st derivative by default:
    k = 1; 
elseif ( k == 0 )
    % Nothing to do here!
    return
end    

% Differentiate with respect to the continuous variable by default:
if ( nargin < 3 )
    dim = 1;
end

%% Take difference across 2nd dimension:
if ( dim == 2 )
    % Differentiate values across dim:
    f.values = diff(f.values, k, dim);
    
    % Differentiate coefficients across dim:
    f.coeffs = diff(f.coeffs, k, dim);
    
    % Tidy up an empty result:
    if ( isempty(f.values) )
        f = f.make(); % Make an empty CHEBTECH.
    end
    
    return    
end

%% Differentiate in dim = 1:

% Get the length:
n = size(f.values,1);

% Get the coefficients:
c = f.coeffs;

% Loop for higher derivatives:
while ( k > 0 )
    % Decrease k
    k = k - 1;
    
    % Derivative of a constant is zero:
    if ( n == 1 )
        f = f.make(zeros(1, size(f.values, 2)), f.vscale, f.epslevel);
        return
    end
    
    % Compute new coefficients using recurrence:
    c = newcoeffs_der(c);
    
    c(abs(c)<f.epslevel) = 0;
    
    % Length of polynomial has decreased by 1:
    n = n - 1;
    
    % Update:
    v = f.chebpolyval(c);
    
% [TODO] If these should be retained, please add a comment explaining. 
%     mu = -1/n*log(f.epslevel);
%     M = f.vscale;
%     f.epslevel =  4*sqrt(sinh(mu)/M)*n*log(n)*f.epslevel*f.vscale;
    
%     rho = exp(mu);
%     phi = sqrt(rho^4+1)/(rho-1).^2;
%     f.epslevel =  (f.epslevel*f.vscale)*(rho-1)*phi*n/sinh(mu)*coth(mu*n)/max(abs(v), [], 1);
    
    f.epslevel = n*log(n)*f.epslevel*f.vscale;
    f.vscale = max(abs(v), [], 1);
end

% Store new coefficients and values:
f.coeffs = c;
f.values = v;

% Compute new values:
% f.values = f.chebpolyval(c);

% [TODO] If these should be retained, please add a comment explaining. 
% % Update the vertical scale:
% % f.vscale = max(f.vscale, max(abs(f.values), [], 1));
% f.vscale = max(abs(f.values), [], 1);

end
      
%%
% Recurrence relation for coefficients of derivative.
function cout = newcoeffs_der(c)
    % C is the vector of coefficients of a Chebyshev polynomial.
    % COUT are the coefficients of its derivative.
% [TODO] These comments should be updated to explain the vectorised case.
    
    [n, m] = size(c);
    cout = zeros(n+1, m);                     % Initialize vector {c_r}
    w = repmat(2*(n-1:-1:1)', 1, m);
    v = [zeros(2, m) ; w.*c(1:end-1,:)];      % Temporal vector
    cout(1:2:end,:) = cumsum(v(1:2:end,:));   % Compute c_{n-2}, c_{n-4},...
    cout(2:2:end,:) = cumsum(v(2:2:end,:));   % Compute c_{n-3}, c_{n-5},...
    cout(end,:) = .5*cout(end,:);             % Adjust the value for c_0
    cout = cout(3:end,:);                     % Trim unneeded coefficients
    
end
