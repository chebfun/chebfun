function h = plus(f, g)
%+   Plus for BALLFUN.
% F + G adds F and G. F and G can be scalars or BALLFUNs.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fIsBallfun = isa(f, 'ballfun');
gIsBallfun = isa(g, 'ballfun');

% Empty check
if isempty( f )
    h = g;
    if ~gIsBallfun
        h = ballfun(h);
    end
    return
end

% Empty check
if isempty( g )
    h = f;
    if ~fIsBallfun
        h = ballfun(h);
    end
    return
end

if (fIsBallfun && gIsBallfun)
    [mf,nf,pf] = size(f); 
    [mg,ng,pg] = size(g); 
    m = max(mf,mg); 
    n = max(nf,ng);
    p = max(pf,pg); 
    X = zeros(m,n,p);
    X(1:mf,floor(n/2)+1-floor(nf/2):floor(n/2)+nf-floor(nf/2),floor(p/2)+1-floor(pf/2):floor(p/2)+pf-floor(pf/2)) = f.coeffs;
    X(1:mg,floor(n/2)+1-floor(ng/2):floor(n/2)+ng-floor(ng/2),floor(p/2)+1-floor(pg/2):floor(p/2)+pg-floor(pg/2)) = ...
    X(1:mg,floor(n/2)+1-floor(ng/2):floor(n/2)+ng-floor(ng/2),floor(p/2)+1-floor(pg/2):floor(p/2)+pg-floor(pg/2)) + g.coeffs;
    h = ballfun(X,'coeffs');
elseif (fIsBallfun && isnumeric(g))
    S = size(f);
    X = f.coeffs;
    % Add the constant g
    X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) = X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) + g;
    h = ballfun(X,'coeffs');
elseif (isnumeric(f) && gIsBallfun)
    S = size(g);
    X = g.coeffs;
    % Add the constant f
    X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) = X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) + f;
    h = ballfun(X,'coeffs');
else
    error('BALLFUN:plus:unknown', ...
          ['Undefined function ''plus'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end
end
