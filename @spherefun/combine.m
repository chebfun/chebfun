function g = combine(g,h)
% COMBINE  Combines two spherefuns together.  
%
% f = combine(g,h) combines g and h into one spherefun.  This is useful
% when g and h have the following properties:
% 
%   g has a CDR decomposition such that C is even and R is pi periodic
%   h has a CDR decomposition such that C is odd and R is pi anti-periodic
%
% If they do not have this property then g+h or plus(g,h) should be used.
%
% See also SPLIT

if ~isa(g,'spherefun') || ~isa(h,'spherefun')
    error('SPHEREFUN:combine:unknown',['Undefined function ''combine'' for ' ...
        'input argument of type %s and %s.'], class(g),class(h));
end

if isempty(g)
    g = h;
    return
elseif isempty(h)
    return
end

idpg = g.idxPlus;
idph = h.idxPlus;
idmg = g.idxMinus;
idmh = h.idxMinus;

pivots = [g.pivotValues(idpg);h.pivotValues(idph);...
          g.pivotValues(idmg);h.pivotValues(idmh)];
cols = [g.cols(:,idpg) h.cols(:,idph) g.cols(:,idmg) h.cols(:,idmh)];
rows = [g.rows(:,idpg) h.rows(:,idph) g.rows(:,idmg) h.rows(:,idmh)];
indices = [g.pivotIndices(:);h.pivotIndices(:)];
locations = [g.pivotLocations;h.pivotLocations];

% Sort the results 
[ignore,perm] = sort( abs(pivots), 1, 'descend' );
pivots = pivots(perm);
cols = cols(:,perm);
rows = rows(:,perm);
indices = indices(perm,:);
locations = locations(perm,:);

% Figure out where the plus/minus terms went
idx = 1:length(pivots);
idx = idx(perm);
plusFlag = idx <= numel([idpg;idph]);
idxPlus = find( plusFlag );
idxMinus = find( ~plusFlag );

% Assemble the results into g.
g.cols = cols;
g.rows = rows;
g.pivotValues = pivots;
g.pivotIndices = indices;
g.pivotLocations = locations;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;
g.nonZeroPoles = g.nonZeroPoles || h.nonZeroPoles;

end

