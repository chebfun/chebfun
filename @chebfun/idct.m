function y = idct(u, type)
%CHEBFUN.IDCT   Inverse discrete cosine transform.
%   CHEBFUN.IDCT(U, TYPE) returns in the inverse discrete cosine transform
%   (inverse DCT) of type KIND on the column vector U. If TYPE is not given it
%   defaults to 3. So far, only types 1, 2, and 3 are supported.
%
%   If U is a matrix, the inverse DCT is applied to each column.
%
% See also CHEBFUN.DCT.

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
        
        % IDCT-I is a scaled DCT-I.
        y = ( 2 / (n-1) ) * chebfun.dct(u, 1);
    
    case 2

        % IDCT-II is a scaled DCT-III.
        y = ( 2 / n ) * chebfun.dct(u, 3); 
    
    case 3
        
        % IDCT-III is a scaled DCT-II:
        y = ( 2 / n ) * chebfun.dct(u, 2);
        
    case 4
        
        % IDCT-IV is a DCT-IV:
        y = ( 2 / n ) * chebfun.dct(u, 4);
        
    otherwise
    
        error('CHEBEFUN:CHEBFUN:idct:type', 'Unknown/unimplemented IDCT type.');
    
end

end