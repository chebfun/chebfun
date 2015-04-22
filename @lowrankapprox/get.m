function val = get( f, propName )
%GET       GET method for LOWRANKAPPROX class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the LOWRANKAPPROX object F. Valid entries for the string PROP are:
%    'DOMAIN'
%    'COLS'
%    'ROWS' 
%    'PIVOTVALUES'
%    'PIVOTLOCATIONS'

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Loop through an array of LOWRANKAPPROX objects.
if ( numel(f) > 1 )
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

% Get the properties.
switch ( propName )
    case 'domain'
        val = f.domain;
    case 'cols'
        val = f.cols;
    case 'rows'
        val = f.rows;
    case 'pivotValues'
        val = f.pivotValues;
    case 'pivotLocations'
        val = f.pivotLocations;
    otherwise
        error('CHEBFUN:LOWRANKAPPROX:get:propName', ...
            [propName,' is not a valid LOWRANKAPPROX property.'])
end
