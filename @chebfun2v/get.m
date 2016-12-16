function val = get( f, propName )
%GET   Get CHEBFUN2V properties.
%
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the CHEBFUN2V F. Valid entries for the string PROP are:
%    'components'  - The CHEBFUN2 in the first component.
%    'nComponents'  - The number of components in a CHEBFUN2.
%    'isTransposed' - Is the CHEBFUN2V a column or row vector?

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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
        error('CHEBFUN:CHEBFUN2V:get:propName', ...
            [propName, ' is not a valid CHEBFUN2V property.'])
end

end
