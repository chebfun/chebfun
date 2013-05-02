function out = sum(g, pref)

% SUM	Definite integral from -1 to 1
% SUM(G) is the integral from -1 to 1 of G.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Unbounded domain map. This works somewhat as domain truncation. For functions
% that decay slowly, this is inaccurate. Exponential decay should work well.

ends = g.domain;
if ( nargin < 2 )
    pref = fun.pref();
end

vends = [get(g, 'lval'), get(g, 'rval')];
tol = 10*get(g, 'epslevel')*get(g, 'vscale');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant case:
if ( length(g) == 1 )
    val = vends(1);
    if ( abs(val) <= 10*pref.fun.eps*vscale );
        out = 0;
    else
        out = inf*sign(val);
    end
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if not zero at infinity: (unbounded integral, simple case)
unbounded = [];
if ( isinf(ends(1)) && abs(vends(1)) > tol )
    unbounded(1) = sign(vends(1))*inf;
end
if ( isinf(ends(2)) && abs(vends(2)) > tol )
    unbounded(2) = sign(vends(2))*inf;
end
if ( ~isempty(unbounded) )
    out = sum(unbounded);
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for sufficient decay at infinity: (1/x is not enough).

% Besides having a zero at (+- 1), the fun should decrease towards the endpoint.
% Decaying faster than 1/x^2 results in a double root, otherwise the integral
% will be unbounded.

% y = g.onefun.points();
% unbounded = [];
% if ( isinf(ends(2)) )
%     gtmp = g.onefun; 
%     gtmp.values = gtmp.values./(1-y);
%     gtmp = extrapolate(gtmp);
%     if ( abs(gtmp(end)) > 1e3*tol && ...
%             diff(gtmp((end-1:end))./diff(y(end-1:end))) > -g.onefun.vscale/g.onefun.hscale )
%         unbounded(1) = inf*sign(g.onefun.values(end-1));
%     elseif ( abs(gtmp(end)) > tol )
%         warning('FUN:sum:slowdecay', ...
%             'Result is likely inaccurate. (Slow decay).')
%     end
% end
% if ( isinf(ends(1)) )
%     gtmp = g.onefun; 
%     gtmp.values = gtmp.values./(1+y);
%     gtmp = extrapolate(gtmp);
%     if ( abs(gtmp(1)) > 1e3*tol && ...
%             diff(gtmp(1:2)./diff(y(1:2))) < g.onefun.vscale/g.onefun.hscale )
%         unbounded(2) = inf*sign(g.onefun.values(2));
%     elseif ( abs(gtmp(1)) > tol )
%         warning('FUN:sum:slowdecay', ...
%             'Result is likely inaccurate. (Slow decay).')
%     end
% end
% if ( ~isempty(unbounded) )
%     out = sum(unbounded);
%     return
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we reach here, the function decays sufficiently fast that we can apply the
% chain rule directly:
% pref.fun.extrapolate = true; 
% pref = onefun.pref(pref, pref.fun);
p = chebtech.pref;
pref.chebtech = p.chebtech;
pref.chebtech.tech = 'cheb1';
g.onefun = onefun.constructor(@(x) feval(g.onefun, x).*g.mapping.der(x), [], [], pref);
out = sum(g.onefun);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


