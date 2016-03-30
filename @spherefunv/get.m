function val = get(f, propName)
%GET   Get SPHEREFUNV properties.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the SPHEREFUNV F. Valid entries for the string PROP are:
%    'components'   - The components of the SPHEREFUN.
%    'isTransposed' - Is the SPHEREFUNV a column or row vector?

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
    case 'isTransposed'
        val = f.isTransposed;
    otherwise
        error('SPHEREFUN:SPHEREFUNV:get:propName', ...
            [propName, ' is not a valid SPHEREFUNV property.'])
end

end
