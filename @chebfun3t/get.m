function val = get(f, propName)
%GET   GET method for CHEBFUN3T class.
%   P = GET(F, PROP) returns the property P specified in the string PROP 
%   from the CHEBFUN F. Valid entries for the string PROP are:
%    'COEFFS'
%    'VSCALE'
%    'DOMAIN'

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Loop through an array of CHEBFUN3 objects.
if ( numel(f) > 1 )
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

% Get the properties.
switch ( propName )
    case 'coeffs'
        val = f.coeffs;
    case 'vscale'
        val = f.vscale;
    case 'domain'
        val = f.domain;
    otherwise
        error('CHEBFUN:CHEBFUN3T:get:propName', ...
            [propName,' is not a valid CHEBFUN3T property.'])
end