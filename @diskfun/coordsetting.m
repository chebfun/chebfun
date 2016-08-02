function iscart = coordsetting( varargin )
% COORDSETTING     Decide if the function is in polar or Cartesian
%   ISCART = COORDSETTING( F ) returns ISCART = 1 if F.COORDS is set to
%   'CART'; otherwise ISCART = 0. 
%  
%   ISCART = COORDSETTING( F, STR ) returns ISCART = 1 if STR = 'cart' and 
%   returns ISCART = 0 if STR = 'polar; 
%
%   ISCART = COORDSETTING( F, ..., STR ) is the same as COORDSETTING( F,
%   STR ).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab diskfun:
f = varargin{ 1 }; 

% First check global setting and apply flag:
if ( strcmpi(f.coords, 'cart') )
    iscart = 1; 
elseif ( strcmpi(f.coords, 'polar') )
    iscart = 0; 
end

% Second search for user-supplied 'polar' flag in arguments: 
isPolar = find(strcmp(varargin, 'polar'));
if ( any( isPolar ) )
    iscart = 0; 
end

% Third search for user-supplied 'cart' setting in feval:
isCart = find(strcmp(varargin, 'cart'));
if ( any( isCart ) )
    iscart = 1; 
end

end




