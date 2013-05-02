function g = cumsum(g)
% CUMSUM	Indefinite integral

dom = g.domain;
vscale = get(g, 'vscale');
vends = [get(g, 'lval'), get(g, 'rval')];
tol = 10*get(g, 'epslevel')*vscale;

% Constant case
if ( length(g) == 1 )
    if ( norm(g, Inf) <= tol )
        g = unbndfun(0, dom);
        return
    else
        error('FUN:cumsum:unbdblow', ...
        ['Representation of functions that blowup on unbounded intervals ' ...
            'has not been implemented in this version'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if not zero at infinity (unbounded integral, simple case)
if ( isinf(dom(1)) && abs(vends(1)) > tol )
    error('FUN:cumsum:unbdblow', ...
        ['Representation of functions that blowup on unbounded intervals ' ...
        'has not been implemented in this version'])
end
if ( isinf(dom(2)) && abs(vends(end)) > tol )
    error('FUN:cumsum:unbdblow', ...
        ['Representation of functions that blowup on unbounded intervals ' ...
        'has not been implemented in this version'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Besides having a zero at (+- 1), the fun should decrease towards the endpoint.
% Decaying faster than 1/x^2 results in a double root.

% if isinf(ends(2))
%     gtmp = g; gtmp.vals = gtmp.vals./(1-y);
%     gtmp = extrapolate(gtmp,pref,y);
%     if abs(gtmp.vals(end)) > 1e3*tol &&  diff(gtmp.vals((end-1:end))./diff(y(end-1:end))) > -g.scl.v/g.scl.h
%         error('FUN:cumsum:unbdblow','Representation of functions that blowup on unbounded intervals has not been implemented in this version')
%     else
%         g.vals(end) = 0;
%         if abs(gtmp.vals(end)) > tol
%             warning('FUN:cumsum:slowdecay','Representation is likely inaccurate')
%         end
%     end
%     
% end
% if isinf(ends(1))
%     gtmp = g; gtmp.vals = gtmp.vals./(1+y);
%     gtmp = extrapolate(gtmp,pref,y);
%     if abs(gtmp.vals(1)) > 1e3*tol && diff(gtmp.vals(1:2)./diff(y(1:2))) < g.scl.v/g.scl.h
%         error('chebfun:cumsum:unbdblow','Representation of functions that blowup on unbounded intervals has not been implemented in this version')
%     else
%         g.vals(1) = 0;
%         if abs(gtmp.vals(1)) > tol
%             warning('FUN:cumsum:slowdecay','Representation is likely inaccurate')
%         end
%     end
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chain rule:
pref.fun.extrapolate = true; 
pref = onefun.pref(pref, pref.fun);
g.onefun = onefun.constructor(@(x) feval(g.onefun, x).*g.mapping.der(x), [], [], pref);
g.onefun = cumsum(g.onefun);

end
