function h = times(f,g)
% .*    Pointwise multiplication for SPHEREFUN objects.
%
%   F.*G multiplies SPHEREFUN objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'spherefun') && isa(g, 'spherefun') )
    % Grady's faster times for rank 1 functions: 
    if ( length( f ) == 1 ) 
        [C, D, R] = cdr( f ); 
        h = g; 
        onesForC = sqrt(abs(D))*ones(1,length(g));
        onesForR = sign(D)*onesForC;
        h.cols = (C*onesForC).*g.cols;
        h.rows = (R*onesForR).*g.rows;
        % Switch the parity terms if needed: 
        % plus*plus = plus -> do nothing
        % plus*minus = minus -> do nothing
        % minus*minus = plus -> switch minus to plus
        % minus*plus = minus -> switch plus to minus
        if ( isempty( f.idxPlus ) )  % f is a minus
            h.idxMinus = g.idxPlus;
            h.idxPlus = g.idxMinus;
        end
        % Both terms have to be non-zero at the poles for the product to be
        % non-zero at the poles.
        h.nonZeroPoles = f.nonZeroPoles & g.nonZeroPoles;
    elseif ( length( g ) == 1 ) 
         h = times(g, f);
    else
        % Let separableApprox deal with the multiplication
        h = times@separableApprox(f,g);
    end
else
    % Let separableApprox deal with the multiplication
    h = times@separableApprox(f,g);
end

end
