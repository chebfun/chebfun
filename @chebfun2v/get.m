function val = get(f,propName)
% GET   Get chebfun2v properties.
%
% P = GET(F,PROP) returns the property P specified in the string PROP from
% the chebfun2v F. Valid entries for the string PROP are:
%   'components'  - The chebfun2 in the first component.
%   'nComponents'  - The number of components in a chebfun2.
%   'isTransposed' - Is the chebfun2v a column or row vector?

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


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
        error('CHEBFUN2v:get:propnam',[propName,' is not a valid chebfun2v property.'])
end

end