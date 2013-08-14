function f = diff(f, k, dim)
%DIFF   Derivative of a CHEBTECH.
%   DIFF(F) is the derivative of F and DIFF(F, K) is the Kth derivative.
%
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference along
%   dimension DIM. For DIM = 1, this is the same as above. For DIM = 2, this
%   is a finite difference along the columns of an array-valued CHEBTECH.
%   If F has L columns, an empty CHEBTECH will be returned for K >= L.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the CHEBTECH G of length n is represented as
%       \sum_{r=0}^{n-1} b_r T_r(x)
% its derivative is represented with a CHEBTECH of length n-1 given by
%       \sum_{r=0}^{n-2} c_r T_r(x)
% where c_0 is determined by
%       c_0 = c_2/2 + b_1;
% and for r > 0,
%       c_r = c_{r+2} + 2*(r+1)*b_{r+1},
% with c_n = c_{n+1} = 0.
%
% [Reference]: Page 34 of Mason & Handscomb, "Chebyshev Polynomials". Chapman &
% Hall/CRC (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check the inputs:

% Trivial case of an empty CHEBTECH:
if ( isempty(f) )
    return
end

if ( nargin < 2 || isempty(k) )
    % Order of derivative not passed in. Assume 1st derivative by default:
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

% If k >= n, we know the result will be the zero function:
if ( k >= n ) 
    z = zeros(size(f, 2));
    f = f.make(z, z, f.hscale);
    return
end

% Loop for higher derivatives:
while ( k > 0 ) % Note that n > k.
    % Decrease k:
    k = k - 1;
    
    % Compute new coefficients using recurrence:
    c = computeDerCoeffs(c);
    
    c(abs(c)<f.epslevel) = 0;
    
    % Length of polynomial has decreased by 1:
    n = n - 1;
    
    % Update:
    v = f.chebpolyval(c);
    
    % Update epslevel and the vertical scale: (See CHEBTECH CLASSDEF file for
    % documentation)
    f.epslevel = n*log(n)*f.epslevel*max(f.vscale); % [TODO]: Vector epslevel?
    f.vscale = max(abs(v), [], 1);
end

% Store new coefficients and values:
f.coeffs = c;
f.values = v;

end
      
function cout = computeDerCoeffs(c)
%COMPUTEDERCOEFFS   Recurrence relation for coefficients of derivative.
%   C is the matrix of Chebyshev coefficients of a (possibly array-valued)
%   CHEBTECH object.  COUT is the matrix of coefficients for a CHEBTECH object
%   whose columns are the derivatives of those of the original.
    
    [n, m] = size(c);
    cout = zeros(n+1, m);                     % Initialize vector {c_r}
    w = repmat(2*(n-1:-1:1)', 1, m);
    v = [zeros(2, m) ; w.*c(1:end-1,:)];      % Temporal vector
    cout(1:2:end,:) = cumsum(v(1:2:end,:));   % Compute c_{n-2}, c_{n-4},...
    cout(2:2:end,:) = cumsum(v(2:2:end,:));   % Compute c_{n-3}, c_{n-5},...
    cout(end,:) = .5*cout(end,:);             % Adjust the value for c_0
    cout = cout(3:end,:);                     % Trim unneeded coefficients
    
end
