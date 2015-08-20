function [f, g, newBreaksLocF, newBreaksLocG] = tweakDomain(f, g, tol, side)
%TWEAKDOMAIN   Adjust nearby common break points in domains of CHEBFUN objects.
%   [F, G] = TWEAKDOMAIN(F, G) adjusts the domains of definition of two CHEBFUN
%   objects F and G if one or more of the entries in F.DOMAIN and G.DOMAIN are
%   sufficiently close. In particular, if abs(F.DOMAIN(j)-G.DOMAIN(k)) <
%   1e-15*max(HSCALE(F),HSCALE(G)) = TOL, then F.DOMAIN(j) and G.DOMAIN(k) are
%   set to (F.DOMAIN(j)+G.DOMAIN(k))/2 (or the nearest integer value if this is
%   less than TOL away). However, if either F or G has two breakpoints which are
%   very close, for example, F.DOMAIN(k+1) - F.DOMAIN(k) < 2*TOL, then these
%   will not be adjusted.
%
%   [F, G] = TWEAKDOMAIN(F, G, TOL) uses the specified tolerance TOL for
%   determining nearby break points.
%
%   [F, G] = TWEAKDOMAIN(F, G, TOL, SIDE) allows control over which of F.domain
%   and G.domain should be given preference for the output.
%       SIDE = 0 : No preference. Take average of F.domain & G.domain (default)
%       SIDE = -1: Give preference to F.domain
%       SIDE = +1: Give preference to G.domain
%
%   [F, G, J, K] = TWEAKDOMAIN(F, G, ...) returns the indices of the modified
%   entries F.DOMAIN(J) and G.DOMAIN(K).
%
%   F = TWEAKDOMAIN(F), where F is a quasimatrix, tweaks the domain of each of
%   the columns of F.
%
%   F = TWEAKDOMAIN(F, DOM), where DOM is numeric or a DOMAIN object tweaks the
%   columns of F against themselves and against DOM. In this case, preference
%   is given to the DOM (i.e., SIDE = 1).
%
% See also CHEBFUN/OVERLAP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse some inputs:
if ( nargin == 1 )
    % Single input is equivalent to TWEAKDOMAIN(F, []);
    g = [];
end
if ( nargin < 3 )
    tol = [];
end
% Return if either f or g are empty as there is nothing to do here:
if ( isempty(f) || isempty(g) )
    newBreaksLocF = [];
    newBreaksLocG = [];
    return
end

if ( isnumeric(g) || isa(g, 'domain' ) )
    % Deal with single input case, where F is typically an quasimatrix. Here we
    % want to ensure that all the columns of F have compatible breaks.
    dom = unique(g);
    if ( numel(f) == 1 && isempty(dom) )
        % Nothing to do in the scalar case if dom is empty.
        return
    end
    if ( numel(dom) >= 2 )   
        domGiven = true;
        dom = g;
    else
        domGiven = false;
        dom = domain(f);                % If not, then use the domain of f.
    end
    if ( isempty(tol) )
    % Set a tolerance relative to the horizontal scale:
        hsg = norm(dom([1, end]), inf);
        % Unbounded domains are defined to have hscale = 1:
        if ( isinf(hsg) ), hsg = 1; end
        tol = 1e-15*max(hscale(f), hsg);
    end
%     dummy = chebfun(0, dom);            % A dummy CHEBFUN to tweak against.
    dummy = struct('domain', dom);      % No - Use struct to avoid overhead!
    newBreaksLocF = cell(1, numel(f));  % Initialise storage.
    for k = 1:numel(f)                  % Loop over columns of f.
        [f(k), dummy, newBreaksLocF{k}, ignored] = ...
            tweakDomain(f(k), dummy, tol, 1);
    end
    if ( domGiven )
        g = dummy.domain;
    end
    return
end

if ( numel(f) > 1 || numel(g) > 1 )
    error('CHEBFUN:CHEBFUN:tweakDomain:quasi', ...
        'tweakDomain does not support multiple quasimatrix inputs.');
end

if ( isempty(tol) )
    % Set a tolerance relative to the horizontal scale:
    hs = max(hscale(f), hscale(g));
    tol = 1e-15*hs;
end

if ( nargin < 4 || isempty(side) )
    side = 0;
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
if ( side == 0 )      % Average
    newBreaks = (f.domain(newBreaksLocF) + g.domain(newBreaksLocG))/2;
elseif ( side < 0 )   % Take f
    newBreaks = f.domain(newBreaksLocF);
else                 % Take g
    newBreaks = g.domain(newBreaksLocG);
end

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
    if ( ~isa(g, 'chebfun') )
        return
    elseif ( k == 1 )                   % First break point is special.
        g.funs{k} = changeMap(g.funs{k}, g.domain([k,k+1]));
    elseif ( k == numel(g.funs) + 1 )   % As is the last.
        g.funs{k-1} = changeMap(g.funs{k-1}, g.domain([k-1,k]));
    else
        g.funs{k-1} = changeMap(g.funs{k-1}, g.domain([k-1,k]));
        g.funs{k} = changeMap(g.funs{k}, g.domain([k,k+1]));
    end
end

end
