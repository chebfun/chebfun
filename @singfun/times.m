function s = times(f,g)
%.* Multiply two singfuns

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Singfuns are superior to doubles and smoothfuns. Those inputs need to be
% promoted to singfuns with zero exponents.
if isa(f,'double')
    %f = smoothfun.constructor(f);
    f = chebtech.constructor(f, vscale, hscale);
elseif isa(g,'double')
    g = smoothfun.constructor(g);
end

if isa(f,'smoothfun')
    f = singfun(f, [0,0]);
elseif isa(g,'smoothfun')
    g = singfun(g, [0,0]);
end

fExps = f.exponents;
gExps = g.exponents;
prefs = singfun.pref('smoothInput',1);

if ( all( abs(fExps-gExps) < tol ) )
    % Case 1: Exponents exactly alike. Just add the smooth parts.
    s = singfun( f.smoothPart + g.smoothPart, f.exponents, prefs );
    
elseif ( abs(round(fExps-gExps) - (fExps-gExps) ) < tol )
    % Case 2: Both exponents differ by integers. Factor out the common
    % singular parts to leave the sum of smooth quotients.
    
    % At each endpoint, the smaller exponent will be factored out
    % of both summands.
    
    % Start off each compensating quotient as 1.
    factorF = @(x) 1;
    factorG = @(x) 1;
    for side = 1:2
        % The smaller of the two exponents is the exponent of the sum.
        [e,k] = sort([ fExps(side),gExps(side)] );
        newExps(side) = e(1);
        
        % The quotient factor is the difference in the exponents.
        if ( side == 1 )
            newFactor = @(x) (1+x).^diff(e);
        else
            newFactor = @(x) (1-x).^diff(e);
        end
        
        % Who had the smaller exponent? The other one gets the factor.
        if ( k(1) == 1 )
            factorG = @(x) factorG(x).*newFactor(x);
        else
            factorF = @(x) factorF(x).*newFactor(x);
        end
    end
    
    % FIXME Do we need to worry about scales here?
    factorF = smoothfun.constructor(factorF);
    factorG = smoothfun.constructor(factorG);
    newSmooth = factorF.*f.smoothPart + factorG.*g.smoothPart;
    s = singfun(newSmooth, newExps, prefs);
    
else
    % Case 3: Nontrivial exponent difference.       
    error('Chebfun:singfun:badAddition',...
        'The resulting sum is not of the form (smooth)x(singular endpoints).')
    
end

end