function y = dst(u, type)
%CHEBFUN.DST   Discrete sine transform.
%   CHEBFUN.DST(U, TYPE) returns in the discrete sine transform (DST) of type
%   KIND on the column vector U. If TYPE is not given it defaults to 1 (which is
%   consistent with MATLAB's PDE toolbox).
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

persistent w1 w2 %#ok<TLEV>

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
        
        % Ensure real/imaginary output for real/imaginary input:
        if ( isreal(u) )
            y = real(y);
        elseif ( isreal(1i*u) )
            y = imag(y);
        end

    case 2
        
        % Equivalent to evaluating a ChebW expansion at 2nd kind points 
        % (up to a diagonal scaling). Also the inverse of DST-III.

        % Weights:  
        if ( size(w1, 1) ~= n+1 )
            w1 = sin( .5*pi*(0:n).'/n );
            w2 = cos( .5*pi*(1:n).'/n );
        end
        W1 = repmat(w1, 1, m);
        W2 = repmat(w2(1:end-1), 1, m);
        z = zeros(1, m);
        
        v1 = chebfun.dct( W1.*[z ; u(1:end-1,:) ; 2*u(end,:)], 1); 
        v2 = chebfun.dst( W2.*u(1:end-1,:), 1);
        y =  [v2 ; z] - v1(2:end,:);
        
    case 3
        
        % Equivalent to evaluating a ChebU expansion at 1st kind points 
        % (up to a diagonal scaling). 

        if ( size(w1, 1) ~= n+1 )
            w1 = sin( .5*pi*(0:n)'/n );
            w2 = cos( .5*pi*(1:n)'/n );
        end
        W1 = repmat(w1, 1, m);
        W2 = repmat(w2, 1, m);
        z = zeros(1, m);

        v1 = W1.*chebfun.dct( [ z ; u ], 1); 
        v2 = W2.*[ chebfun.dst( u(1:end-1,:), 1) ; z ];
        y = v2 - v1(2:end,:);
        
    case 4 
        % Equivalent to evaluating a ChebW expansion at 1st kind points 
        % (up to a diagonal scaling).

        y = chebfun.dct( u(end:-1:1,:) , 4);
        y(2:2:end,:) = -y(2:2:end,:);
                
    case 5
        % Relate DSTV of length N to a DSTI of length 2N: 
        y = chebfun.dst( [ u ; zeros(n, m)], 1);
        y = y(2:2:end,:); 
        
    case 6
        % Relate DSTVI of length N to a DSTII of length 2N-1: 
        y = zeros(2*n, m); 
        y(1:2:end, :) = u; 
        y = chebfun.dst( y, 1);
        y = y(1:n,:); 
        
    case 7
        % Relate DSTVIII of length N to a DSTIII of length 2N-1: 
        y = chebfun.dst( [ u ; zeros(n, m)], 1);
        y = y(1:2:end,:); 
        
    case 8
        % Relate DSTVIII of length N to a DSTIV of length 2N-1: 
        y = chebfun.dst( [ u ; zeros(n+1, m)], 3);
        y = y(1:2:2*n,:); 
        
    otherwise
    
        error('CHEBEFUN:CHEBFUN:dct:type', 'Unknown DST type.');
        
end

end
