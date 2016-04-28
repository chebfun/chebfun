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
            [mf, nf] = length(f); 
            [mg, ng] = length(g);
            m = max(mf, mg);
            n = max(nf, ng); 
            [Cf, Df, Rf] = coeffs2(f, m, n); 
            [Cg, Dg, Rg] = coeffs2(g, m, n);
            Cf = chebtech2.coeffs2vals( Cf ); 
            Rf = chebtech2.coeffs2vals( Rf ); 
            Cg = chebtech2.coeffs2vals(Cg)*sqrt(Dg);
            Rg = chebtech2.coeffs2vals(Rg)*sqrt(Dg); 
            H = zeros(m, n); 
            for k = 1:size(Cg,2)
               H = H + (bsxfun(@times, Cf, Cg(:,k)))*Df*(bsxfun(@times, Rf, Rg(:,k)))'; 
            end
            h = chebfun2( H );
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
