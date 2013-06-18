function f = compose(f, op, g, pref)
%COMPOSE   Compostition of BNDFUN objects.
%   H = COMPOSE(F, OP) returns a BNDFUN representing OP(F) where F is also a
%   BNDFUN object, and OP is a function handle.
%
%   H = COMPOSE(F, OP, G) returns OP(F, G) where F and G are BNDFUN objects, and
%   OP is a function handle.
%
%   H = COMPOSE(F, G) returns a BNDFUN representing G(F), where both F and G are
%   also BNDFUN objects. If the range of F is not enclosed by the domain of G,
%   an error is thrown. Notice that the domain of H will be the same as the
%   domain of F.
%
%   H = COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), or COMPOSE(F, G, PREF)
%   uses the options passed by the prefences structure PREF. In particular, one
%   can pass a PREF.(class(F)).refinmentFunction which takes advantage of the
%   fact that F (and possibly OP or G) are BNDFUN objects.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin == 2 )
    if ( isa(op, 'bndfun') )
        % op is a bndfun! Rename it to g for clarity
        g = op;
        
        % For composition to work with two BNDFUN objects, the range of F must
        % be within the domain of G. Check whether that is the case.
        
        % [TODO]: Should BNDFUN compose assume that it's being called with
        % sensible inputs? Presumably, we're having to match up pieces at a
        % higher level, and multiple calls to find mins and maxes will be
        % expensive.
        
        % The tolerance of the match is determined by the hscale of F and G.
        tol = max(get(f,'hscale'), get(g,'hscale'))*eps;
        
        % min and max values of F
        mvals = minandmax(f);
        fMinVal = min(mvals(1,:));
        fMaxVal = max(mvals(2,:));
        
        % Left and right endpoints of the domain of G
        dom = g.domain;
        gA = dom(1);
        gB = dom(2);
        
        if ( (fMinVal < (gA - tol)) || (fMaxVal > (gB + tol)) )
            % Throw an error if domains don't match:
            error('BNDFUN:compose:domainMismatch', ...
                ['For composition G(F), the range of F must be enclosed ', ...
                'by the domain of G.'])
        end
        
        % Map the range of F to a range which lies in [-1 1]. Note that we have
        % already verified that the domain of G, that is [gA,gB], encloses the
        % range of F, that is [fMinVals, fMaxVals]. Hence, this guarantees that
        % the values of fMapped, that is, its range, will be within [-1,1],
        % which is necessary, as the domain of the onefun of G is [-1,1].
        %
        % However, notice that it is not necessary that the range of fMapped is
        % the entire interval [-1,1], it only has to be included in [-1,1]. This
        % mapping is a linear transformation of the range of the ONEFUN of f
        % (which in fact is the same as the range of F) to a (subset) of the
        % domain of the ONEFUN of G, i.e., [-1,1].
        fMapped = (2*f.onefun - (gA + gB))./(gB - gA);
        
        % Do the composition
        f.onefun = compose(fMapped, g.onefun);    
    else
        % OP is an operator!
        f.onefun = compose(f.onefun, op);
    end
    
elseif ( nargin == 3 )
    
    if ( isstruct(g) )
        % Third argument passed was a preference structure.
        f.onefun = compose(f.onefun,op, [], g);
    else
        % Third argument passed was a bndfun. Compose the onefun of f with the
        % onefun of g:
        f.onefun = compose(f.onefun,op,g.onefun);
    end
else
    f.onefun = compose(f.onefun,op,g.onefun,pref);
end

end
