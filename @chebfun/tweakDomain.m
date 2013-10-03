function [f, g, newBreaksLocF, newBreaksLocG] = tweakDomain(f, g, tol)
%TWEAKDOMAIN   Adjust nearby common break points in domains of CHEBFUN objects.
%   [F, G] = TWEAKDOMAIN(F, G) adjusts the domains of definition of two CHEBFUN
%   objects F and G if one or more of the entries in F.DOMAIN and G.DOMAIN are
%   sufficiently close. In particular, if 
%     abs(F.DOMAIN(j)-G.DOMAIN(k)) < 1e-15*max(HSCALE(F),HSCALE(G)) = TOL,
%   then F.DOMAIN(j) and G.DOMAIN(k) are set to (F.DOMAIN(j)+G.DOMAIN(k))/2 (or
%   the nearest integer value if this is less than TOL away). However, if either
%   F or G has two breakpoints which are very close, for example, F.DOMAIN(k+1)
%   - F.DOMAIN(k) < 2*TOL, then these will not be adjusted.
%
%   [F, G] = TWEAKDOMAIN(F, G, TOL) uses the specified tolerance TOL for
%   determining nearby break points.
%
%   [F, G, J, K] = TWEAKDOMAIN(F, G) returns the indices of the modified
%   entries F.DOMAIN(J) and G.DOMAIN(K).
%
% See also CHEBFUN/OVERLAP.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Return if either f or g are empty as there is nothing to do here:
if ( isempty(f) || isempty(g) )
    newBreaksLocF = [];
    newBreaksLocG = [];
    return
end

if ( nargin < 3 )
    % Set a tolerance relative to the horizontal scale:
    hs = max(hscale(f), hscale(g));
    tol = 1e-15*hs;
end

%% Decide which breaks need adjusting:
fDom = f.domain;
gDom = g.domain;

% We don't want to change breakpoints in F or G that are too close:
tinyIntsF = diff(fDom) < 2*tol; 
fDom([tinyIntsF, 0] | [0, tinyIntsF]) = NaN; % Set to NaN to fool abs(a-b) in
tinyIntsG = diff(gDom) < 2*tol;              % diffDomains below. This prevents
gDom([tinyIntsG, 0] | [0, tinyIntsG]) = NaN; % these points from being included.

% Find where there are two break points which are _almost_ the same
diffDomains = bsxfun(@(a,b) abs(a - b), fDom, gDom.');
idx = diffDomains > 0 & diffDomains < tol;
newBreaksLocF = max(idx, [], 1);
newBreaksLocG = max(idx, [], 2).';

%% Determine new breaks:

% Take the average of the nearby break points to create a new one:
newBreaks = (f.domain(newBreaksLocF) + g.domain(newBreaksLocG))/2;
% Move this to an integer if it is sufficiently close:
rndIdx = abs(round(newBreaks) - newBreaks) < tol;
newBreaks(rndIdx) = round(newBreaks(rndIdx));

%% Update F and G:

% Update the .domain entries in the CHEBFUN objects:
f.domain(newBreaksLocF) = newBreaks;
g.domain(newBreaksLocG) = newBreaks;

% Change the maps used by FUN objects on either side of new breakpoints:
% NOTE: This is slightly inefficient if consecutive break points are being
% adjusted, as the map for one of the FUN objects is being changed twice.
% However, this is a fairly inexpensive operation, so we don't worry too much.
newBreaksLocF = find(newBreaksLocF);
newBreaksLocG = find(newBreaksLocG);
for k = newBreaksLocF
    % Change the FUN objects belonging to f:
    if ( k == 1 )                       % First break point is special.
        f.funs{k} = changeMap(f.funs{k}, f.domain([k,k+1]));
    elseif ( k == numel(f.funs) + 1 )   % As is the last.
        f.funs{k-1} = changeMap(f.funs{k-1}, f.domain([k-1,k]));
    else
        f.funs{k-1} = changeMap(f.funs{k-1}, f.domain([k-1,k]));
        f.funs{k} = changeMap(f.funs{k}, f.domain([k,k+1]));
    end
end
for k = newBreaksLocG
    % Change the FUN objects belonging to g:
    if ( k == 1 )                       % First break point is special.
        g.funs{k} = changeMap(g.funs{k}, g.domain([k,k+1]));
    elseif ( k == numel(g.funs) + 1 )   % As is the last.
        g.funs{k-1} = changeMap(g.funs{k-1}, g.domain([k-1,k]));
    else
        g.funs{k-1} = changeMap(g.funs{k-1}, g.domain([k-1,k]));
        g.funs{k} = changeMap(g.funs{k}, g.domain([k,k+1]));
    end
end

end
