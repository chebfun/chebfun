function val = get( f, propName )
%GET GET method for BALLFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the BALLFUN F. Valid entries for the string PROP are:
%    'coeffs', 'isReal', 'domain'

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the properties.
switch ( propName )
    case 'coeffs'
        val = f.coeffs;
    case 'isReal'
        val = f.isReal;
    case 'domain'
        val = f.domain;     
otherwise
        error('CHEBFUN:BALLFUN:get:propName', ...
            [propName,' is not a valid BALLFUN property.'])
end