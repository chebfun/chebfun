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

n = size(u, 1);

switch type
    
    case 1
        
        % IDCT-I is a scaled DCT-I.
    
        u = flipud(u);
        y = 2/(n-1)*chebfun.dct(u, 1);
        y = flipud(y);
    
    case 2
        
        % IDCT-II is a scaled DCT-III.
        
        u = sqrt(2/n)*u;
        u(1,:) = sqrt(2)*u(1,:);
        y = chebfun.dct(u, 3); 
        y = flipud(y);
    
    case 3
        
        % IDCT-III is a scaled DCT-II:

        u = flipud(u);
        y = chebfun.dct(u, 2);
        y = sqrt(2/n)*y;
        y(1,:) = sqrt(2)*y(1,:);
        
    otherwise
    
        error('CHEBEFUN:CHEBFUN:idct:type', 'Unknown/unimplemented IDCT type.');
    
end

end