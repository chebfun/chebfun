function y = idct(u, type)
%CHEBFUN.IDCT   Inverse discrete cosine transform.
%   CHEBFUN.IDCT(U, TYPE) returns in the inverse discrete cosine transform
%   (inverse DCT) of type KIND on the column vector U. If TYPE is not given it
%   defaults to 2. So far, types 1-4 are supported.
%
%   If U is a matrix, the inverse DCT is applied to each column.
%
%   IDCTs are scaled in many different ways. We have decided to be consistent
%   with Wikipedia: http://en.wikipedia.org/wiki/Discrete_cosine_transform.
%
%   Note that the above means that CHEBFUN.IDCT(U) is not the same as IDCT(U),
%   where IDCT(U) is the implementation in the Matlab signal processing toolbox.
%   The two are related by 
%       IDCT(U) = CHEBFUN.DCT(E*U)
%   where E = sqrt(2)*eye(n); E(1,1) = 2;
%
% See also CHEBFUN.DCT, CHEBFUN.DST, CHEBFUN.IDST.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Implement IDCT5-8. 

% Default to kind 3:
if ( nargin < 2 )
    type = 2;
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