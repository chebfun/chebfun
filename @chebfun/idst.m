function y = idst(u, type)
%CHEBFUN.IDST   Inverse discrete sine transform.
%   CHEBFUN.IDST(U, TYPE) returns in the inverse discrete sine transform
%   (inverse DST) of type KIND on the column vector U. If TYPE is not given it
%   defaults to 1 (which is consistent with Matlab's PDE toolbox). So far, only
%   types 1-4 are supported.
%
%   If U is a matrix, the inverse DST is applied to each column.
%
%   IDSTs are scaled in many different ways. We have decided to be consistent
%   with Wikipedia: http://en.wikipedia.org/wiki/Discrete_sine_transform.
%
%   Note that CHEBFUN.IDST(U) is the same as IDST(U) from Matlab's PDE toolbox.
%
% See also CHEBFUN.DST, CHEBFUN.DCT, CHEBFUN.IDCT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Implement IDST5-8. 

% Default to kind 1:
if ( nargin < 2 )
    type = 1;
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
    
        error('CHEBEFUN:CHEBFUN:idct:type', 'Unknown/unimplemented IDST type.');
    
end

end