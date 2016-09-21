function val = get(f, propName)
%GET   Get CHEBFUN3V properties.
%   P = GET(F, PROP) returns the property P specified in the string PROP 
%   from the CHEBFUN3V object F. Valid entries for the string PROP are:
%   'components'   - CHEBFUN3 components of F.
%   'nComponents'  - The number of components in a CHEBFUN3.
%   'isTransposed' - Is the CHEBFUN3V a column or row vector?

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

switch propName
    case 'components'
        val = f.components;
    case 'nComponents'
        val = f.nComponents; 
    case 'isTransposed'
        val = f.isTransposed;
    otherwise
        error('CHEBFUN:CHEBFUN3V:get:propName', ...
            [propName, ' is not a valid CHEBFUN3V property.'])
end

end