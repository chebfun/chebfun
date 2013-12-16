function A = subsasgn(A, sr, B)
% Replace blocks in a chebmatrix.
% A(...) = B replaces the block or submatrix referenced on the left side of
% the equality with the blocks given on the right side.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if isa(B,'chebmatrix')
    data = B.blocks;
else
    data = { B }; 
end

switch(sr(1).type)
    case {'()', '{}'}
        A.blocks = subsasgn( A.blocks, sr, data);
    otherwise
        A.(sr(1).subs) = B;
end

end