function y = dst(u, type)
%CHEBFUN.DST   Discrete sine transform.
%   CHEBFUN.DST(U, TYPE) returns in the discrete sine transform (DST) of type
%   KIND on the column vector U. If TYPE is not given it defaults to 2. So far,
%   only types 1-4 are supported.
%
%   If U is a matrix, the DST is applied to each column.
%
% See also CHEBFUN.IDST.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Check normalization and document in help text.

% Default to kind 2.
if ( nargin < 2 )
    type = 2;
end

[n, m] = size(u);

switch type
    
    case 1
        
        % Equivalent to evaluating a ChebU expansion at 2nd kind points 
        % (up to a diagonal scaling).
        Z = zeros(1, m); 
        y = fft( [ Z ; u ; Z ; -u(end:-1:1,:)] ); 
        y = imag( y( 2*n+2:-1:n+3, : )/2 );

    case 2
        
        % Equivalent to evaluating a ChebW expansion at 2nd kind points 
        % (up to a diagonal scaling). Also the inverse of DST-III.  
        y = u;
        for k = n - 1 : -1 : 1
            y(k, :) = y(k, :) - y(k+1, :); 
        end
        y = chebfun.dst(y, 4);
        y = bsxfun( @times, y, 2*cos(pi/2/n*((0:n-1)'+1/2)) );
        
%         y = chebfun.dct( u(end:-1:1,:) , 2);
%         y = bsxfun(@times, y, (-1).^(0:n-1)' );
        
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
    
        error('CHEBEFUN:CHEBFUN:dct:type', 'Unknown/unimplemented DCT type.');
        
end

end
