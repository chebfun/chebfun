function [funs, ends] = constructor(op, domain, pref)
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.
 
% Set hscale and vscale:
hscale = norm(domain,inf);
if ( isinf(hscale) )
    hscale = 1; 
end
vscale = 0;

% Initial setup.
ends = domain;
nfuns = numel(ends)-1;

% ----------------------------SPLITTING OFF-------------------------------

% In off mode, seek only one piece with length no greater than maxdegree (default is 2^16)
if ( ~pref.chebfun.splitting )
     maxn = pref.chebfun.maxdegree + 1;
     funs{nfuns} = fun.constructor();
     for k = 1:nfuns
         endsk = ends(k:k+1);
         if ( iscell(op) )
             opk = op{k};
         else
             opk = op; 
         end
         [funs{k}, ishappy, vscale] = getfun(opk, endsk, hscale, vscale, pref);
         if ( ~ishappy )
            warning('CHEBFUN:auto', ['Function not resolved, using %d pts.', ...
                ' Have you tried ''splitting on''?'], maxn);
         end
         
     end
     return
end

% ----------------------------SPLITTING ON--------------------------------

pref.chebfun.extrapolate = true;
pref.chebfun.maxSamples = pref.chebfun.splitdegree;
% We extrapolate when splitting so that we can construct functions like
% chebfun(@sign,[-1 1]), which otherwise would not be happy at x = 0.

funs{nfuns} = fun.constructor();
ishappy = ones(1,numel(ends)-1);
% Try to get one smooth piece for the entire interval before splitting interval
for k = 1:nfuns
     [funs{k}, ishappy(k), vscale] = getfun(op, ends(k:k+1), hscale, vscale, pref);
end
sad = ~ishappy;

% MAIN LOOP If the above didn't work, enter main loop and start splitting (at
% least one breakpoint will be introduced).

while any(sad)       
    % If a fun is sad in a subinterval, split this subinterval.
    j = find(sad, 1, 'first');
    a = ends(j);
    b = ends(j+1);
    
    [edge, vscale] = chebfun.detectedge(op, [a b], hscale, vscale, pref);
    
%     htol = 1e-14*hscale;
%     if ( abs(a - edge) <= htol )
%         edge = a + diff([a b])/100;
%     elseif ( abs(b - edge) <= htol )
%         edge = b - diff([a b])/100;
%     end
    
    % Try to obtain happy child funs on each new subinterval.
    % ------------------------------------

    %  left
    [child1, happy1, vscale] = getfun(op, [a, edge], hscale, vscale, pref);

    %  right
    [child2, happy2, vscale] = getfun(op, [edge, b], hscale, vscale, pref);
    
    % Insert to existing funs:
    funs = [funs(1:j-1), {child1, child2}, funs(j+1:end)];
    ends = [ends(1:j), edge, ends(j+1:end)];

    % Check for happiness:
    sad = [sad(1:j-1), ~happy1, ~happy2, sad(j+1:end)];
    if ( numel(ends) > 40 )
%         error('too many intervals!')
        warning('too many intervals!')
        return
    end
   
end

end

function [g, ishappy, vscale] = getfun(op, domain, hscale, vscale, pref)
%GETFUN controls the construction of funs

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Initial setup
htol = 1e-14*hscale;
% If the interval is very small skip adaptation and return a constant
% (This should never be happen, though!)
if ( diff(domain) < 2*htol )
    [g, ishappy, vscale] = smallInterval(op, domain, hscale, vscale, pref);
    return
end

% Call the fun constructor.
pref = fun.pref(pref, pref.chebfun);
g = fun.constructor(op, domain, hscale, vscale, pref);
ishappy = get(g, 'ishappy');

% Update the vertical scale.
vscale = max(vscale, get(g, 'vscale'));

end

function [g, ishappy, vscale] = smallInterval(op, domain, hscale, vscale, pref)
    pref = fun.pref(pref, pref.chebfun);
    g = fun.constructor(op(mean(domain)), domain, hscale, vscale, pref);
    ishappy = get(g, 'ishappy');
    vscale = max(vscale, get(g,'vscale'));
end


    