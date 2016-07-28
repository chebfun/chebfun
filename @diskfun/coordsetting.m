function iscart = coordsetting(varargin)
% COORDSETTING: Given input to a function, determine whether the
% function evaluation should occur with respect to polar or Cartesian coordinates. 
% The first input element must be a diskfun. 

%figure out if Cartesian or polar
f = varargin{1}; 

%first check global setting and apply
if strcmpi(f.coords,'cart')
    iscart =1; 
elseif strcmpi(f.coords, 'polar')
    iscart =0; 
end

%second search for user-supplied 'polar' setting in arguments: 
isPolar = find(strcmp(varargin,'polar'));
if ( any(isPolar) )
    iscart = 0; 
end
%third search for user-supplied 'cart' setting in feval:
isCart = find(strcmp(varargin,'cart'));
if ( any(isCart) )
    iscart = 1; 
end

end




