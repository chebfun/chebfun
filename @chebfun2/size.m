function varargout = size( f , dim)
% Size   Size of a chebfun2. 

if ( ( nargin == 1 ) && ( nargout <= 1 ) ) 
        varargout = {[Inf, Inf]}; 
elseif ( ( nargout == 2 ) && ( nargout <=1 ) )
        varargout = {Inf, Inf}; 
elseif ( nargin == 2 ) && ( ( dim == 1 ) || ( dim == 2 ) )
    varargout = {Inf}; 
else
    error('CHEBFUN2:SIZE:OUTPUTS','Too many output arguments.');
end

end