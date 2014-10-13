function y = dct(u, type)
%CHEBFUN.DCT   Discrete cosine transform.
%   CHEBFUN.DCT(U, TYPE) returns in the discrete cosine transform (DCT) of type
%   KIND on the column vector U. If TYPE is not given it defaults to 2. So far,
%   types 1-4 are supported.
%
%   If U is a matrix, the DCT is applied to each column.
%
%   DCTs are scaled in many different ways. We have decided to be consistent
%   with Wikipedia: http://en.wikipedia.org/wiki/Discrete_cosine_transform.
%
%   Note that the above means that CHEBFUN.DCT(R) is not the same as DCT(U),
%   where DCT(U) is the implementation in the Matlab signal processing toolbox.
%   The two are related by 
%       DCT(U) = E*CHEBFUN.DCT(U)
%   where n = size(U, 1) and E = sqrt(2/n)*speye(n); E(1,1) = 1/sqrt(n).
%
% See also CHEBFUN.IDCT, CHEBFUN.DST, CHEBFUN.IDST.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Implement DCT5-8. 

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
        
        y = fft( [ u ; u(end-1:-1:2 ,: ) ]/2 ); 
        y = y( 1:n, : );
        
    case 2
        
        % Equivalent to evaluating a ChebV expansion at 2nd kind points 
        % (up to a diagonal scaling). Also the inverse of DCT-III 
        % (up to a diagonal scaling) and this is how it is implemented.  
        
        u = flipud( u );
        y = ( n / 2 ) * chebtech1.vals2coeffs( u );        
        y(1,:) = 2 * y(1, :); 
        
    case 3
        
        % Equivalent to evaluating a ChebT expansion at 1st kind points 
        % (up to a diagonal scaling). Implemented using the connection to a 
        % real-even DFT of half-shifted output, see CHEBTECH2.COEFFS2VALS().  
        
        u(1,:) = .5*u(1, :); 
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
