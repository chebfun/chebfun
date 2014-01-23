function A = createBlocks(disc,item)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin < 2 )
    item = disc.source;
end
if isa(item,'linBlock')
    disc.source = item;
    A = disc.source.stack( disc );
elseif isa(item,'chebfun')
    A = disc.toValues(item);
elseif isnumeric(item)
    A = item;
else
    error('Unrecognized block type.')
end
end
