function varargout = quiver3( F, varargin )
%QUIVER3   Quiver plot of a SPHEREFUNV.
%   Same as QUIVER.
%
% See also QUIVER.

h = quiver(F,varargin{:});

if ( nargout > 0 )
    varargout = {h};
end

end
