function A = subsasgn(A, sr, B)
% Change properties in a chebmatrix.
% A(...) = B replaces the block or submatrix referenced on the left side of
% the equality with the blocks given on the right side.
%
% A.(...) = B changes a named property of A. To change a preference value, use a
% double dot reference. For example,
%
%  A.prefs.discretization = @ultraS;
%  A.prefs.maxTotalLength = 2000;

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

switch(sr(1).type)
    
    case {'()', '{}'}
        if isa(B,'chebmatrix')
            data = B.blocks;
        else
            data = { B };
        end

        A.blocks = subsasgn( A.blocks, sr, data);
    
    otherwise  % dot reference (property)
        
        if ( length(sr) > 1 )
            % Nested ref. This occurs when setting a pref, most notably.
            p = A.(sr(1).subs);
            A.(sr(1).subs) = subsasgn( p, sr(2:end), B);
        else
            A.(sr(1).subs) = B;
        end

end

end