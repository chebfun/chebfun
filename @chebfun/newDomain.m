function g = newDomain(g, newDom)
%NEWDOMAIN   Change of domain of a CHEBFUN.
%  NEWDOMAIN(G, DOM) returns the CHEBFUN G but moved to the domain DOM. This is
%  done with a linear map. DOM may be a vector of length G.ends, or a two-vector
%  (in which case all breakpoints are scaled by the same amount).

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Unbounded domains.

% Current breakpoints:
oldDom = g.domain;

if ( numel(newDom) == numel(oldDom) )
    % All new breakpoints are given!
elseif ( numel(newDom) == 2 )
    % Scale breakpoints:
    c = oldDom(1); 
    d = oldDom(end);
    a = newDom(1);  
    b = newDom(2);
    newDom = (b - a)*(oldDom - c)/(d - c) + a;
else
    error('CHEBFUN:CHEBFUN:newDomain:numints', 'Inconsistent domains.');
end

for k = 1:numel(g.funs)
    % Update the domains of each of the funs:
    g.funs{k} = changeMap(g.funs{k}, newDom(k:k+1));
end

% Update the CHEBFUN:
g.domain = newDom;

end
