function y = idst(u, type)
%CHEBFUN.IDST   Inverse discrete sine transform.
%   CHEBFUN.IDST(U, TYPE) returns in the inverse discrete sine transform
%   (inverse DST) of type KIND on the column vector U. If TYPE is not given it
%   defaults to 1 (which is consistent with Matlab's PDE toolbox). So far
%   types 1-7 are supported.
%
%   If U is a matrix, the inverse DST is applied to each column.
%
%   IDSTs are scaled in many different ways. We have decided to be consistent
%   with Wikipedia: http://en.wikipedia.org/wiki/Discrete_sine_transform.
%
%   Note that CHEBFUN.IDST(U) is the same as IDST(U) from Matlab's PDE toolbox.
%
% See also CHEBFUN.DST, CHEBFUN.DCT, CHEBFUN.IDCT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to kind 1:
if ( nargin < 2 )
    type = 1;
end

[n, m] = size(u); %#ok<ASGLU>

switch type
    
    case 1
        
        % IDST-I is a scaled DST-I.
        y = ( 2 / (n+1) ) * chebfun.dst(u, 1);
    
    case 2

        % IDST-II is a scaled DST-III.
        y = ( 2 / n ) * chebfun.dst(u, 3); 
    
    case 3
        
        % IDST-III is a scaled DST-II:
        y = ( 2 / n ) * chebfun.dst(u, 2);
        
    case 4
        
        % IDST-IV is a DST-IV:
        y = ( 2 / n ) * chebfun.dst(u, 4);
        
    case 5
        
        % IDST-V is a DST-V:
        y = ( 4 / (2*n+1) ) * chebfun.dst(u, 5);
        
    case 6
        
        % IDST-VI is a DST-VII:
        y = ( 4 / (2*n+1) ) * chebfun.dst(u, 7);
        
    case 7
        
        % IDST-VII is a DST-VI:
        y = ( 4 / (2*n+1) ) * chebfun.dst(u, 6);
        
    case 8
        
        error('CHEBFUN:CHEBFUN:IDST:EIGHT', 'Not implemented')
%         y = NaN;
        
    otherwise
    
        error('CHEBFUN:CHEBFUN:idct:type', 'Unknown IDST type.');
    
end

end
