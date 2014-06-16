function varargout = chebkind(kind)
%CHEBKIND   Set the default Chebyshev grid type.
%   CHEBKIND 1 or CHEBKIND 2 globally sets the default grid used by the CHERBFUN
%   constructor to be Chebyshev points of the first or second kind,
%   respectively, and is equivalent to setting the default 'tech' preference in
%   CHEBFUNPREF to be either @chebtech1 or @chebtech2.
%
%   CHEBKIND by itself returns the current default Chebyshev grid type. If the
%   default technology is not Chebyshev, then zero is returned.
%
% See also CHEBTECH1, CHEBTECH2, CHEBFUNPREF, CHEBFUNPREF.SETDEFAULTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )
    tech = feval(chebfunpref().tech);
    if ( isa(tech, 'chebtech') )
        tech = class(tech);
        kind = str2double(tech(end));
    else
        kind = 0;
    end
elseif ( (kind == 1) || any(strcmp(kind, {'1', '1st'})) )
    chebfunpref.setDefaults('tech', @chebtech1);
    kind = 1;
elseif ( (kind == 2) || any(strcmp(kind, {'2', '2nd'})) )
    chebfunpref.setDefaults('tech', @chebtech2);
    kind = 2;
else
    error('CHEBFUN:chebkind:unknown', ...
        'Unknown input to CHEBKIND(). Valid options are 1 or 2.');
end

if ( nargout > 0 )
    varargout{1} = kind;
end

end