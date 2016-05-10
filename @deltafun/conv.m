function h = conv(f, g)
%CONV   Convolution of DELTAFUN objects.
%   H = CONV(F, G) produces the convolution of DELTAFUN objects F and G:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [a + c, b + d]
%                    /
%   The output H is a cell array of DELTAFUNs or CLASSICFUNS. If the result of
%   the convolution of F and G is zero, an empty cell is returned. The cell
%   array returned is used by higher level convolutions.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = deltafun();
    return
end

if ( ~isa(f, 'deltafun') )
    deltaMagF = [];
    deltaLocF = [];
    funF = f;
else
    f = simplifyDeltas(f);
    if ( isa(f, 'deltafun') )
        deltaMagF = f.deltaMag;
        deltaLocF = f.deltaLoc; 
        funF = f.funPart;
    else
        deltaMagF = [];
        deltaLocF = [];
        funF = f;
    end
end

if ( ~isa(g, 'deltafun') )
    deltaMagG = [];
    deltaLocG = [];
    funG = g;
else
    g = simplifyDeltas(g);
    if ( isa(g, 'deltafun') )
        deltaMagG = g.deltaMag;
        deltaLocG = g.deltaLoc; 
        funG = g.funPart;
    else
        deltaMagG = [];
        deltaLocG = [];
        funG = g;
    end
end

% Extract the domains of f and g:
domF = funF.domain;
domG = funG.domain;
a = domF(1);
b = domF(end);
c = domG(1);
d = domG(end);

% Get the threshold for deltafunctions:
pref = chebfunpref();
deltaTol = pref.deltaPrefs.deltaTol;

% Compute the convolution of funParts and append it to the output cell:
h = conv(funF, funG);

% Remove zero funs for simplicity:
zeroIndices = cellfun( @(hk) iszero(hk), h);
h(zeroIndices) = [];

%% Get all the deltafunction contributions
% (f + df) * (g + dg) = f*g + dg * (f + df) + df * (g + dg) - df * dg
%                     = f*g + dg * (f + df/2) + df * (g + dg/2)

% Contributions due to deltafunctions in F:
% df * (g + dg/2):
h = convDeltas(deltaMagF, deltaLocF, g, c, d, deltaTol, h);

% Contributions due to delta functions in G:
% dg * (f + df/2);
h = convDeltas(deltaMagG, deltaLocG, f, a, b, deltaTol, h);

% Check for emptiness:
if ( isempty(h) )
    h = {fun.constructor(0, struct('domain', [a, b]))};
    return
end

%% Make sure that h is a cell array with non-overlapping funs.

% First extract the breakpints from funs forming h:
doms = cellfun(@(hk) hk.domain, h, 'uniformoutput', 0);
breakPoints = sort(unique(horzcat(doms{:})));

% Coalesce breakpoints that are extremely close together.
if ( any(isinf(breakPoints)) )
    tol = 2*eps;  % hscale of functions on unbounded domains is 1.
else
    tol = 2*eps*norm(breakPoints, Inf);
end

% TODO:  Take the average of points which are close together instead of just
% picking one more or less arbitrarily?
closeBreaks = abs(diff(breakPoints)) < tol;
breakPoints(closeBreaks) = [];

% Restrict each FUN to lie within the new breakpoints.
nFuns = length(h);
H = {};
deltaTol = pref.deltaPrefs.proximityTol;
for k = 1:nFuns
    hk = h{k};
    dom = hk.domain;

    % The domain endpoints of each computed FUN must lie _exactly_ within the
    % set of breakpoints or RESTRICT will fail.
    idx1 = find(abs(breakPoints - dom(1)) < tol);
    idx2 = find(abs(breakPoints - dom(2)) < tol);
    hk = changeMap(hk, breakPoints([idx1, idx2]));

    % If we just changed the map of a DELTAFUN, the locations of deltas at the
    % domain endpoints also need to be updated.
    if ( isa(hk, 'deltafun') )
        deltaLoc = hk.deltaLoc;
        domk = hk.domain;
        if( ~isempty(deltaLoc) )
            if ( abs(deltaLoc(1) - domk(1) ) < deltaTol )
                deltaLoc(1) = domk(1);
            end
            if ( abs(deltaLoc(end) - domk(2)) < deltaTol )
                deltaLoc(end) = domk(2);
            end
            hk.deltaLoc = deltaLoc;
        end
    end

    % Do the restriction.
    hk = restrict(hk, breakPoints(idx1:idx2));
    if ( ~isa(hk, 'cell') )
        hk = {hk};
    end
    H = [H, hk]; %#ok<AGROW>
end

% If H has the same number of funs, nothing further needs to be done:
if ( length(H) == length(breakPoints) - 1 )
    return
end

% H has more funs as a result of restrict. Re-initialize h with 0 BNDFUNS:
nFuns = length(breakPoints) - 1;
h = {};
for k = 1:nFuns
    data.domain = breakPoints(k:k+1);
    h = [h, {bndfun(0, data)}]; %#ok<AGROW>
end

% Loop through H and add each fun to the relevant domain piece in h:
for k = 1:length(H)
    hk = H{k};
    dom = hk.domain;
    if ( abs(dom(2) - dom(1)) > tol )
        idx = find(breakPoints == dom(1));
        h{idx} = h{idx} + hk; %#ok<AGROW>
    end
end

end

function h = convDeltas(deltaMagF, deltaLocF, g, c, d, deltaTol, h)
%   H = CONVDELTAS() convolves the function G with deltafunctions described by
%   deltaMagF and deltaLocF. G originally lives on [c, d] and the output is
%   returned in H. This function is called twice so deltafunctions in g are 
%   scaled by half to get the correct delta-delta contributions.

%%
% Contributions due to deltafunctions in F:
% df * (g + dg/2):
[m, n] = size(deltaMagF);
% Loop through the delta function matrix:
for i = 1:m
    for j = 1:n
        if ( abs(deltaMagF(i, j)) > deltaTol )
            % Take appropriate derivative, scale and shift the function:
            hij = deltaMagF(i, j) * changeMap(diff(g,i-1), deltaLocF(j) + [c d]);
            % The half below is to make sure that delta-delta interaction
            % is not counted twice:
            if ( isa(hij, 'deltafun') )
                hij.deltaMag = hij.deltaMag/2;                
            end
            
            % If the result is non-zero, append it:
            if ( ~iszero(hij) )
                h = [h, {hij}]; %#ok<AGROW>
            end
        end
    end
end
end
