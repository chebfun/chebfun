function y = idst(u, type)
%CHEBFUN.IDST   Inverse discrete sine transform.
%   CHEBFUN.IDST(U, TYPE) returns in the inverse discrete sine transform
%   (inverse DST) of type KIND on the column vector U. If TYPE is not given it
%   defaults to 3. So far, only types 1, 2, and 3 are supported.
%
%   If U is a matrix, the inverse DST is applied to each column.
%
% See also CHEBFUN.DST.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Check normalization and document in help text.

% Default to kind 3:
if ( nargin < 2 )
    type = 3;
end

[n, m] = size(u);

switch type
    
    case 1
        
        % IDST-I is a scaled DST-I.
        y = ( 2 / (n+1) ) * chebfun.dst(u, 1);
    
    case 2

        % IDST-II is a scaled DST-III.
        u(n, :) = 2 * u(n, :);
        y = ( 2 / n ) * chebfun.dst(u, 3); 
        y(n, :) = y(n, :) / 2;
    
    case 3
        
        % IDST-III is a scaled DST-II:
        u(n, :) = u(n, :) / 2;
        y = ( 2 / n ) * chebfun.dst(u, 2);
        y(n, :) = 2 * y(n, :);
        
    case 4
        
        % IDST-IV is a DST-IV:
        y = ( 2 / n ) * chebfun.dst(u, 4);
        
    otherwise
    
        error('CHEBEFUN:CHEBFUN:idct:type', 'Unknown/unimplemented IDCT type.');
    
end

end