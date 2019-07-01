function val = get( F, propName )
%GET       GET method for BALLFUNV class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the BALLFUNV object F. Valid entries for the string PROP are:
%    'comp'
%    'components'

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the properties.
switch ( propName )
    case 'comp'
        val = F.comp;
    case 'components'
        val = F.comp;
otherwise
        error('CHEBFUN:BALLFUNV:get:propName', ...
            [propName,' is not a valid BALLFUNV property.'])
end