function y = dct(u, type)
%CHEBFUN.DCT   Discrete cosine transform.
%   CHEBFUN.DCT(U, TYPE) returns in the discrete cosine transform (DCT) of type
%   KIND on the column vector U. If TYPE is not given it defaults to 2. So far,
%   only types 1, 2, and 3 are supported.
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

n = length(u);

switch type
    
    case 1
        
        u([1,end],:) = .5*u([1,end],:);
        u = flipud(u);
        y = chebtech2.coeffs2vals(u);
    
    case 2
        
        u = flipud(u);
        y = sqrt(n/2)*chebtech1.vals2coeffs(u);
        y = flipud(y);
        y(1,:) = sqrt(2)*y(1,:);
        
    case 3
    
        u(1,:) = .5*u(1,:);
        u = flipud(u);
        y = chebtech1.coeffs2vals(u);    
    
    otherwise
    
        error('CHEBEFUN:CHEBFUN:dct:type', 'Unknown/unimplemented DCT type.');
    
end

end
