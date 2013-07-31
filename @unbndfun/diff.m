function f = diff(f, k, dim)
%DIFF	Derivative of a BNDFUN.
%   DIFF(F) is the derivative of F and DIFF(F, K) is the K-th derivative.
%
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference along
%   dimension DIM. For DIM = 1, this is the same as above. For DIM = 2, this is
%   a finite difference along the columns of an array-valued UNBNDFUN. If F has
%   L columns, an UNBNDFUN with an empty ONEFUN property will be returned for K
%   >= L.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% First derivative by default:
if ( nargin < 2 )
    k = 1;
end

if ( k == 0 )
    % Nothing to do here!
    return
end

if ( nargin == 3 && dim == 2 )
    
    % Difference across columns:
    f.onefun = diff(f.onefun, k, 2);
    
else
    % Differentiate along columns.
    
    % Create an UNBNDFUN of the reciprocal of the derivative of the map:
    gprcp = f; 
    gprcp.onefun = onefun.constructor(@(x) 1./f.mapping.der(x));

    % Loop for higher derivatives:
    for j = 1:k

        % Differentiate the ONEFUN:
        f.onefun = diff(f.onefun);

        % Apply the chain rule:
        f = f.*gprcp;

    end

end

end