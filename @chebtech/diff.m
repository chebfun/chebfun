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

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

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

if ( (nargin < 2) || isempty(k) )
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

if ( dim == 1 )
    % Take difference across 1st dimension:
    f = diffContinuousDim(f, k);
else
    % Take difference across 2nd dimension:
    f = diffFiniteDim(f, k);
end

end

function f = diffFiniteDim(f, k)
% Take kth difference across 2nd dimension (i.e., across columns).

    if ( k >= size(f, 2) )
        % The output will be an empty CHEBTECH:
        f = f.make();
    else 
        for j = 1:k
            % Differentiate coefficients across columns:
            f.coeffs = diff(f.coeffs, 1, 2);
            % Update vscale and epslevel as in PLUS():
            ev = f.epslevel.*f.vscale;
            for l = 1:(size(f, 2)-1)
                f.epslevel(l) = ev(l) + ev(l+1);
            end
            % We've lost a column, so we lose an epslevel:
            f.epslevel(end) = [];
            % New vscale and epslevel:
            vscale = getvscl(f);
            f.epslevel = f.epslevel./vscale;
            f.vscale = vscale;
        end
    end
end

function f = diffContinuousDim(f, k)
% Differentiate in the first dimension (i.e., df/dx).
    
    % Get the coefficients:
    c = f.coeffs;

    % Get their length:
    n = size(c, 1);

    % If k >= n, we know the result will be the zero function:
    if ( k >= n ) 
        z = zeros(size(f, 2));

        data.vscale = z;
        data.hscale = f.hscale;
        f = f.make(z, data);
        return
    end
    
    % Loop for higher derivatives:
    for m = 1:k
        % Compute new coefficients using recurrence:
        c = computeDerCoeffs(c);
        n = n - 1;
    end
    
    % Store new coefficients:
    f.coeffs = c;
    
    % Update epslevel and the vertical scale: (See CHEBTECH CLASSDEF file
    % for documentation)
    newVscale = getvscl(f);
    epslevelBnd = (n*log(n)).^k*(f.epslevel.*f.vscale)./newVscale;
    f.epslevel = updateEpslevel(f, epslevelBnd);
    f.vscale = newVscale;
    
end
      
function cout = computeDerCoeffs(c)
%COMPUTEDERCOEFFS   Recurrence relation for coefficients of derivative.
%   C is the matrix of Chebyshev coefficients of a (possibly array-valued)
%   CHEBTECH object.  COUT is the matrix of coefficients for a CHEBTECH object
%   whose columns are the derivatives of those of the original.

    [n, m] = size(c);
    cout = zeros(n-1, m);                       % Initialize vector {c_r}
    w = repmat(2*(1:n-1)', 1, m);
    v = w.*c(2:end,:);                          % Temporal vector
    cout(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:));   % Compute c_{n-2}, c_{n-4},...
    cout(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:));   % Compute c_{n-3}, c_{n-5},...
    cout(1,:) = .5*cout(1,:);                   % Adjust the value for c_0
end
