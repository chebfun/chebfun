function y = idct(u, type)
%CHEBFUN.IDCT   Inverse discrete cosine transform.
%   CHEBFUN.IDCT(U, TYPE) returns in the inverse discrete cosine transform
%   (inverse DCT) of type KIND on the column vector U. If TYPE is not given it
%   defaults to 2.
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

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Default to kind 2:
if ( nargin < 2 )
    type = 2;
end

[n, m] = size(u); %#ok<ASGLU>

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
        
    case 5 
        
        % IDCT-V is a DCT-V: 
        u(1,:) = u(1,:)/2;
        y = ( 4 / (2*n-1) ) * chebfun.dct(u, 5);
        y(1,:) = y(1,:)/2; 
        
    case 6
        
        % IDCT-VI is a DCT-VII: 
        u(1,:) = u(1,:)/2; 
        y = ( 4 / (2*n-1) ) * chebfun.dct(u, 7);
        y(end,:) = y(end,:)/2;
        
    case 7
        
        % IDCT-VII is a DCT-VI: 
        u(end,:) = u(end,:)/2; 
        y = ( 4 / (2*n-1) ) * chebfun.dct(u, 6);
        y(1,:) = y(1,:)/2; 
        
    case 8 
        
        % IDCT-VIII is a DCT-VIII: 
        y = ( 4 / (2*n+1) ) * chebfun.dct(u, 8);
        
    otherwise
    
        error('CHEBEFUN:CHEBFUN:idct:type', 'Unknown IDCT type.');
    
end

end
