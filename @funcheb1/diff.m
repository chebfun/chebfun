function f = diff(f, k, dim)
%DIFF	Differentiation of a FUNCHEB1
%   DIFF(F) is the derivative of F and DIFF(F, K) is the K-th derivative.
%
%   See also SUM, CUMSUM.

% [TODO]:  Document the DIM argument.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the fun G of length n is represented as
%       SUM_{r=0}^{n-1} C_r T_r(x)
% its derivative is represented by a fun of length n-1 given by
%       SUM_{r=0}^{n-2} c_r T_r (x)
% where c_0 is determined by
%       c_0 = c_2/2 + C_1;
% and for r > 0,
%       c_r = c_{r+2} + 2*(r+1)*C_{r+1},
% with c_{n} = c_{n+1} = 0.
% (See "Chebyshev Polynomials" by Mason and Handscomb, CRC 2002, p. 34.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Check the inputs

% Trivial case of an empty fun:
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

%% Take difference across 2nd dimension (i.e., dim = 2):
if ( dim == 2 )
    % Differentiate values across dim:
    f.values = diff(f.values, k, dim);
    
    % Differentiate coeffs across dim:
    f.coeffs = diff(f.coeffs, k, dim);
    
    % Tidy up an empty result:
    if ( isempty(f.values) )
        f = funcheb1();
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
    
    % Derivative of a constant is zero
    if ( n == 1 )
        f.values = 0;
        f.vscale = 0; 
        f.coeffs = 0;
        return
    end
    
    % Compute new coefficients using recurrence:
    c = newcoeffs_der(c);
    
    % Length of polynomial has decreased by 1:
    n = n-1;
end

% Compute new values:
f.values = funcheb1.chebpolyval(c);

% Store new coefficients:
f.coeffs = c;

% Update the vertical scale:
f.vscale = max(f.vscale, max(f.values(:)));

end
      
%%
% Recurrence relation for coefficients of derivative
function cout = newcoeffs_der(c)
    % C is the coefficients of a Chebyshev polynomial.
    % COUT are the coefficients of its derivative.
    
    [n, m] = size(c);
    cout = zeros(n+1,m);                     % Initialize vector {c_r}
    w = repmat(2*(n-1:-1:1)', 1, m);
    v = [zeros(2,m) ; w.*c(1:end-1,:)];      % Temporal vector
    cout(1:2:end,:) = cumsum(v(1:2:end,:));  % Compute c_{n-2}, c_{n-4},...
    cout(2:2:end,:) = cumsum(v(2:2:end,:));  % Compute c_{n-3}, c_{n-5},...
    cout(end,:) = .5*cout(end,:);            % Rectify the value for c_0
    cout = cout(3:end,:);                    % Trim unneeded coeffs
    
end

