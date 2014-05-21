function A = subsasgn(A, sa, B)
%SUBSASGN  Change a property of a chebmatrix.
%   A(...) = B replaces the block or submatrix referenced on the left side
%   of the equality with the blocks given on the right side.
%
%   A.(...) = B changes a named property of A. 

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

switch(sa(1).type)
    
    case {'()', '{}'}
        if isa(B, 'chebmatrix')
            data = B.blocks;
        else
            data = { B };
        end

        A.blocks = subsasgn( A.blocks, sa, data);
    
    otherwise  % dot reference (property)
        
        if ( length(sa) > 1 )
            % Nested ref. This occurs when setting a pref, most notably.
            p = A.(sa(1).subs);
            A.(sa(1).subs) = subsasgn( p, sa(2:end), B);
        else
            A.(sa(1).subs) = B;
        end

end

end
