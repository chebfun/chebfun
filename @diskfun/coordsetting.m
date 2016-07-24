function iscart = coordsetting(varargin)
% Given input to a function, determines whether the coordinate setting 
% is polar or
% cartesian.
% first element input must be a diskfun. 

%figure out if cartesian or polar
f = varargin{1}; 
%f = f{1};
%vars = varargin{:};
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




