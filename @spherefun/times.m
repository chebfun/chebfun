function h = times(f, g)
% .*   SPHEREFUN multiplication.
%   F.*G multiplies SPHEREFUN objects F and G. Alternatively F or G could be a
%   double.

if ( isa(f, 'spherefun') )    % SPHEREFUN .* ???
    
    if ( isa(g, 'double') )  % SPHEREFUN .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'spherefun') )
        if ( bol )
            h = spherefun(@(th, lam) feval(f, th, lam).*feval(g, th, lam) );
        else
            error('SPHEREFUN:times:domain', 'Inconsistent domains');
        end
    else
        error('SPHEREFUN:times:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end
    
else
    h = times(g, f);
end

end
