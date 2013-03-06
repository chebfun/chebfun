function [edge, vscale] = detectedge(op, domain, hscale, vscale, pref)
    % Parse the inputs:
    if ( nargin < 5 )
        pref = chebfun.pref;
    end
    if ( nargin < 2 )
        domain = pref.chebfun.domain;
    end
    if ( nargin < 3 )
        hscale = norm(domain,inf);
        if ( isinf(hscale) )
            hscale = 1; 
        end
    end
    if ( nargin < 4 )
        vscale = 0;
    end

    % Call detectedge at the fun level:
    [edge, vscale] = fun.detectedge(op, domain, hscale, vscale);
end