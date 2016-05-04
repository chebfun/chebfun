function h = times(f, g)
% .*    Pointwise multiplication for SEPARABLEAPPROX objects.
%
%   F.*G multiplies SEPARABLEAPPROX objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
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
                h = g; 
                onesForC = sqrt(abs(D))*ones(1,length(g));
                onesForR = sign(D)*onesForC;
                h.cols = (C*onesForC).*g.cols;
                h.rows = (R*onesForR).*g.rows;
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
