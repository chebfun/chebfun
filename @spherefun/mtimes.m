function h = mtimes(f, g)
%*	SPHEREFUN multiplication.
%   c*F or F*c multiplies a SPHEREFUN F by a scalar c.
%
% See also TIMES.

if ( isa(f, 'spherefun') )           % SPHEREFUN * ???
    
    if ( isa(g, 'double') )         % SPHEREFUN * DOUBLE
        if ( numel(g) == 1 )
            h = f;
            h.BlockDiag = h.BlockDiag .* g;
        else
            error('SPHEREFUN:mtimes:size', 'Sizes are inconsistent.');
        end
               
    else
        error('SPHEREFUN:mtimes:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
        
    end
    
else
    
    h = mtimes(g.', f.').';
    
end

end
