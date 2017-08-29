function disc = extractBlock(disc, j, k)
%EXTRACTBLOCK   Extract the j-k block from a discretization.
%   DISCJK = EXTRACTBLOCK(DISC, J, K) extracts information relating to the
%   J-K block of the dscretization DISC and returns a new discretization
%   DISCJK. DISCJK.dimension will be modified by the amount specified in
%   DISC.dimAdjust(j,k).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract the j-k block:
if ( iscell(disc.source.blocks) )
    disc.source.blocks = disc.source.blocks{j,k};
end

% Extract the dimension adjustment for this block and adjust the dimension:
if ( numel(disc.dimAdjust) > 1 )
    disc.dimAdjust = disc.dimAdjust(j,k);
end
disc.dimension = disc.dimension + disc.dimAdjust;

% Set the dimension adjustment to zero (as dimension was adjusted above)
disc.dimAdjust = 0;

% Extract the projection Order:
if ( numel(disc.projOrder) > 1 )
    disc.projOrder = disc.projOrder(j);
end

end
