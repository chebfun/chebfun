function y = dct(u, type)
%CHEBFUN.DCT   Discrete cosine transform.
%   CHEBFUN.DCT(U, TYPE) returns in the discrete cosine transform (DCT) of type
%   KIND on the column vector U. If TYPE is not given it defaults to 2. So far,
%   only types 1-4 are supported.
%
%   If U is a matrix, the DCT is applied to each column.
%
% See also CHEBFUN.IDCT.

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
        
        % Equivalent to evaluating a ChebT expansion at 2nd kind points 
        % (up to a diagonal scaling). Implemented using the connection to 
        % a real-even DFT, see CHEBTECH2.COEFFS2VALS().
        
        u([1,end], :) = .5 * u([1,end], :);
        u = flipud( u );
        y = chebtech2.coeffs2vals( u );
        y = flipud( y ); 
        
    case 2
        
        % Equivalent to evaluating a ChebV expansion at 2nd kind points 
        % (up to a diagonal scaling). Also the inverse of DCT-III 
        % (up to a diagonal scaling) and this is how it is implemented.  
        
        u = flipud( u );
        y = ( n / 2 ) * chebtech1.vals2coeffs( u );
        y = flipud( y );
        y(1,:) = 2 * y(1, :); 
        
    case 3
        
        % Equivalent to evaluating a ChebT expansion at 1st kind points 
        % (up to a diagonal scaling). Implemented using the connection to a 
        % real-even DFT of half-shifted output, see CHEBTECH2.COEFFS2VALS().  
        
        u(1,:) = .5*u(1, :);
        u = flipud( u );
        y = chebtech1.coeffs2vals( u );    
        y = flipud( y ); 
        
    case 4 
        
        % Equivalent to evaluating a ChebV expansion at 1st kind points 
        % (up to a diagonal scaling).
        
        y = bsxfun( @times, u, cos(pi/2/n*((0:(n-1))'+1/2)) );
        y = chebfun.dct( y, 2 );
        y(2:n, :) = 2 * y(2:n, :);
        for k = 2 : n
            y(k, :) = y(k, :) - y(k-1, :); 
        end

    otherwise
    
        error('CHEBEFUN:CHEBFUN:dct:type', 'Unknown/unimplemented DCT type.');
        
end

end
