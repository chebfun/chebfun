function varargout = show(c)
% SHOW Display a chebgui in the GUI window.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/chebfun/ for Chebfun information.

[h1, h2] = chebguiWindow(c);

if ( nargout > 0 )
    varargout{1} = h1;
end

if ( nargout > 1 )
    varargout{2} = h2;
end

end

