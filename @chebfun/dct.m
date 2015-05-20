function y = dct(u, type)
%CHEBFUN.DCT   Discrete cosine transform.
%   CHEBFUN.DCT(U, TYPE) returns in the discrete cosine transform (DCT) of type
%   KIND on the column vector U. If TYPE is not given it defaults to 2.
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

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to kind 2.
if ( nargin < 2 )
    type = 2;
end

[n, m] = size(u);

switch type
    
    case 1
        
        % Equivalent to evaluating a ChebT expansion at 2nd kind points 
        % (up to a diagonal scaling). Implemented using the connection to 
        % CHEBTECH2.COEFFS2VALS().

        u([1, end],:) = .5*u([1, end],:);
        y = chebtech2.coeffs2vals(u);
        y = y(end:-1:1,:);
        
    case 2
        
        % Equivalent to evaluating a ChebV expansion at 2nd kind points 
        % (up to a diagonal scaling). Also the inverse of DCT-III 
        % (up to a diagonal scaling) and this is how it is implemented.  
        
        y = ( n / 2 ) * chebtech1.vals2coeffs( u(end:-1:1,:) );        
        y(1,:) = 2 * y(1,:); 
        
    case 3
        
        % Equivalent to evaluating a ChebT expansion at 1st kind points 
        % (up to a diagonal scaling). Implemented using the connection to a 
        % real-even DFT of half-shifted output, see CHEBTECH1.COEFFS2VALS().  
        
        u(1,:) = .5*u(1,:); 
        y = chebtech1.coeffs2vals( u );    
        y = y(end:-1:1,:); 
        
    case 4 
        
        % Equivalent to evaluating a ChebV expansion at 1st kind points 
        % (up to a diagonal scaling).

        v = zeros(2*n, m); 
        v(2:2:end, :) = u; 
        y = chebfun.dct( v, 3 );
        y = y(1:n, :);

    case 5 
        
        % Relate DCTV of length N to a DCTI of length 2N: 
        y = chebfun.dct( [2*u(1, :) ; u(2:end, :) ; zeros(n, m)], 1 ); 
        y = y(1:2:end,:);
        
    case 6 
        
        % Relate DCTVI of length N to a DCTII of length 2N-1: 
        y = chebfun.dct( [u ; zeros(n-1, m)], 2 ); 
        y = y(1:2:end,:);
        
    case 7 
        
        % Relate DCTVII of length N to a DCTIII of length 2N-1: 
        v = zeros(2*n-1, m); 
        v(1:2:end,:) = u; 
        v(1,:) = 2*v(1,:);
        y = chebfun.dct( v, 3 );
        y = y(1:n,:);
        
    case 8 
            
        % Relate DCTVIII of length N to a DCTII of length 2N+1:
        y = chebfun.dct( [u ; zeros(n+1,m)], 2 );
        y = y(2:2:end,:);
        
    otherwise
    
        error('CHEBEFUN:CHEBFUN:dct:type', 'Unknown DCT type.');
        
end

end
