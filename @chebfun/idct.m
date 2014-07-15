function y = idct(u, type)
%CHEBFUN.IDCT   Inverse discrete cosine transform.
%   CHEBFUN.IDCT(U, TYPE) returns in the inverse discrete cosine transform
%   (inverse DCT) of type KIND on the column vector U. If TYPE is not given it
%   defaults to 3.
%
%   If U is a matrix, the inverse DCT is applied to each column.
%
% See also CHEBFUN.DCT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to kind 3:
if ( nargin < 2 )
    type = 3;
end

switch type
    
    case 1
    
        u = flipud(u);
        n = size(u, 1);
        y = 2/(n-1)*chebfun.dct(u, 1);
        y = flipud(y);
    
    case 2
        
        % TODO: Use DCT-III
        y = flipud(chebtech1.vals2coeffs(u));
        y(1,:) = 2*y(1,:);
    
    case 3
        
        error('IDCT Type-III not implemented.');
        
    otherwise
    
        error('unknown kind');
    
end

end