function out = innerProduct(f, g)
%INNERPRODUCT Compute the inner product of two DELTAFUN objects.
%   INNERPRODUCT(F, G) is the inner-product of DELTAFUN F and G. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% Trivial cases:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

% Make sure both arguments are DELTAFUNS:
if ( ~isa(g, 'deltafun') )
    g = deltafun(f, [], []);
end

if ( ~isa(f, 'deltafun') )
    f = deltafun(f, [], []);
end

%%
% The innerProduct has four contributions:
% < (f + df) , (g + dg ) > = <f, g> + <df, g> + <dg, f> + <df, dg>

% <dg, dg>: If two delta functions are at the same location, compute the 
% appropriate infinity and return:
[fcomIdx, gcomIdx] = sameDeltaLocs(f, g);
if ( any(fcomIdx) )
    % If there is an overlap of locations, extract the columns corresponding to
    % the first overlap:
    fIdx = find(fcomIdx, 1);
    gIdx = find(gcomIdx, 1);
    fd = f.deltaMag(:, fIdx);
    gd = g.deltaMag(:, gIdx);
    
    % Extract, the first non-zero product term:
    pref = chebfunpref();
    tol = pref.deltaPrefs.deltaTol;
    fd = fd(find(abs(fd) > tol, 1));
    gd = gd(find(abs(gd) > tol, 1));
    % Assign a signed infinity:
    out = inf * sign(fd*gd);
    return;
end

% Delta functions don't overlap: compute the components of the inner product:
funIP = innerProduct(f.funPart, g.funPart);           % <fFun,gFun>
dfIP = deltaInnerProduct(g, f.deltaMag, f.deltaLoc);  % <g, df>
dgIP = deltaInnerProduct(f, g.deltaMag, g.deltaLoc);  % <f, dg>
out = funIP + dfIP + dgIP;
end


function deltaIP = deltaInnerProduct(g, deltaMag, deltaLoc)    
%DELTAINNERPRODUCT Compute the inner product of a DELTAFUN G with dleta
%functions.

% Handle the empty case:
if ( isempty(deltaLoc) )
    deltaIP = 0;
    return
end

% Compute the derivatives needed:
m = size(deltaMag, 1);
maxDiffOrder = m-1;
G = zeros(m, length(deltaLoc));
G(1,:) = feval(g, deltaLoc);
for k = 1:maxDiffOrder
    g = diff(g);
    G(k+1, :) = feval(g, deltaLoc);
end

% The output is always a scalar double:
deltaIP = 0;
v = ones(maxDiffOrder+1,1); v(2:2:end) = -1;
for k = 1:length(deltaLoc)
    % Apply the definition of inner product with delta functions:
    % <sum(ai dirac^(i)(x-xk)), f(x)> = sum ai*(-1)^i f^(i)(xk)
    ipk = (v.*deltaMag(:, k)).' * G(:, k) ;
    deltaIP = deltaIP + ipk;
end

end