function h = times(f, g)
% .*    Pointwise multiplication for SEPARABLEAPPROX objects.
%
%   F.*G multiplies SEPARABLEAPPROX objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'separableApprox') )    % SEPARABLEAPPROX .* ???
    
    if ( isa(g, 'double') )  % SEPARABLEAPPROX .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'separableApprox') )
        bol = domainCheck(f, g);
        if ( bol )
            % Grady's faster times for rank 1 functions: 
            if ( length( f ) == 1 ) 
                [C, D, R] = cdr( f ); 
                C = C*sqrt(D); 
                R = R*sqrt(D);
                h = g; 
                for k = 1:length( g ) 
                    h.cols(:,k) = C.*(g.cols(:,k)); 
                    h.rows(:,k) = R.*(g.rows(:,k)); 
                end

            elseif ( length( g ) == 1 ) 
                 h = times(g, f);
            else
                % Give up, call the constructor: 
                h = compose(f, @times, g); 
            end
        else
            error('CHEBFUN:SEPARABLEAPPROX:times:domain', 'Inconsistent domains');
        end
    else
        error('CHEBFUN:SEPARABLEAPPROX:times:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end
    
else
    h = times(g, f);
end

end
