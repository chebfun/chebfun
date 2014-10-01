function y = dst(u, type)
%CHEBFUN.DST   Discrete sine transform.
%   CHEBFUN.DST(U, TYPE) returns in the discrete sine transform (DST) of type
%   KIND on the column vector U. If TYPE is not given it defaults to 1 (which is
%   consistent with Matlab's PDE toolbox). So far, only types 1-4 are supported.
%
%   If U is a matrix, the DST is applied to each column.
%
%   DSTs are scaled in many different ways. We have decided to be consistent
%   with Wikipedia: http://en.wikipedia.org/wiki/Discrete_sine_transform.
%
%   Note that CHEBFUN.DST(U) is the same as DST(U) from Matlab's PDE toolbox.
%
% See also CHEBFUN.IDST, CHEBFUN.DCT, CHEBFUN.IDCT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Implement DST5-8. 

% Default to kind 1:
if ( nargin < 2 )
    type = 1;
end

[n, m] = size(u);

switch type
    
    case 1
        
        % Equivalent to evaluating a ChebU expansion at 2nd kind points 
        % (up to a diagonal scaling).
        Z = zeros(1, m); 
        y = fft( [ Z ; u ; Z ; -u(end:-1:1,:)] ); 
        y = -1i*( y( 2*n+2:-1:n+3, : )/2 );

    case 2
        
        % Equivalent to evaluating a ChebW expansion at 2nd kind points 
        % (up to a diagonal scaling). Also the inverse of DST-III.  
        y = u;
        for k = n - 1 : -1 : 1
            y(k, :) = y(k, :) - y(k+1, :); 
        end
        y = chebfun.dst(y, 4);
        y = bsxfun( @times, y, 2*cos(pi/2/n*((0:n-1)'+1/2)) );
        
    case 3
        
        % Equivalent to evaluating a ChebU expansion at 1st kind points 
        % (up to a diagonal scaling).  
        u(n, :) = .5 * u(n, :);
        y = bsxfun(@times, u, 2*cos(pi/2/n*((0:n-1)'+1/2)) );
        y = chebfun.dst( y, 4 );
        for k = 2 : n
            y(k, :) = y(k, :) - y(k-1, :); 
        end
        
    case 4 
        % Equivalent to evaluating a ChebW expansion at 1st kind points 
        % (up to a diagonal scaling).

        y = chebfun.dct( u(end:-1:1,:) , 4);
        y = bsxfun(@times, y, (-1).^(0:n-1)' );
        
    otherwise
    
        error('CHEBEFUN:CHEBFUN:dct:type', 'Unknown/unimplemented DST type.');
        
end

end
