function g = combine(g, h)
%COMBINE    Combines two SPHEREFUN objects together.
%
%   f = combine(g, h) combines g and h into one SPHEREFUN, where g and h
%   have the following properties:
%   g has a CDR decomposition such that C is even and R is pi periodic,
%   h has a CDR decomposition such that C is odd and R is pi anti-periodic.
%
%   If they do not have this property then g+h or plus(g, h) should be used.
%
% See also PARTITION

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( ~isa(g,'spherefun') || ~isa(h,'spherefun') )
    error('CHEBFUN:SPHEREFUN:combine:unknown',['Undefined function ''combine'' '...
        'for input argument of type %s and %s.'], class(g), class(h));
end

if ( isempty(g) )
    g = h;
    return
elseif ( isempty(h) )
    return
end

idPlusG = g.idxPlus;
idPlusH = h.idxPlus;
idMinusG = g.idxMinus;
idMinusH = h.idxMinus;

% Only combine spherefuns that have one strict type of parity.
if ( ~isempty(idPlusG) && ~isempty(idMinusG)) || (~isempty(idPlusH) && ...
        ~isempty(idMinusH) )
    error('CHEBFUN:SPHEREFUN:combine:parity',['Inputs must have oposite ' ...
        'parity. Consider using plus']);
end

pivots = [ g.pivotValues(idPlusG); h.pivotValues(idPlusH); ...
           g.pivotValues(idMinusG); h.pivotValues(idMinusH) ];
cols = [ g.cols(:, idPlusG) h.cols(:, idPlusH) g.cols(:, idMinusG) ...
    h.cols(:, idMinusH) ];
rows = [ g.rows(:, idPlusG) h.rows(:, idPlusH) g.rows(:, idMinusG) ...
    h.rows(:, idMinusH) ];
locations = [ g.pivotLocations; h.pivotLocations ];

numPlus = length(idPlusG) + length(idPlusH);
idxPlus = 1:numPlus;
numMinus = length(idMinusG) + length(idMinusH);
idxMinus = (numPlus+1):(numPlus+numMinus);
 
% TODO: The code below sorts columns and rows according to the size of the
% pivots. This may mess up the assumed location of the non-zero pole term
% for subsequent operations. Therefore we are commenting it out.  We may
% decide later that sorting is a good idea, in which case we will have to
% figure out how to deal with the non-zero pole.
% % Sort the results 
% [ignore, perm] = sort(abs(pivots), 1, 'descend');
% pivots = pivots(perm);
% cols = cols(:, perm);
% rows = rows(:, perm);
% indices = indices(perm, :);
% locations = locations(perm, :);
% 
% % Figure out where the plus/minus terms went
% idx = 1:length(pivots);
% idx = idx(perm);
% plusFlag = idx <= numel([idPlusG; idPlusH]);
% idxPlus = find(plusFlag);
% idxMinus = find(~plusFlag);

% Assemble the results into g:
g.cols = cols;
g.rows = rows;
g.pivotValues = pivots;
g.pivotLocations = locations;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;
g.nonZeroPoles = g.nonZeroPoles || h.nonZeroPoles;

end