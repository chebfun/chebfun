function g = combine(g,h)
%COMBINE   Combines two diskfuns together.  
%
%   f = combine(g,h) combines g and h into one diskfun.  This is useful
%   when g and h have the following properties:
% 
%     g has a CDR decomposition such that C is even and R is pi periodic
%     h has a CDR decomposition such that C is odd and R is pi anti-periodic
%
%   If they do not have this property then g+h or plus(g,h) should be used.
%
% See also DISKFUN/PARTITION.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ~isa(g,'diskfun') || ~isa(h,'diskfun')
    error('CHEBFUN:DISKFUN:combine:unknown',['Undefined function ''combine'' for ' ...
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

% Only combine diskfuns that have one strict type of parity.
if (~isempty(idpg) && ~isempty(idmg)) || (~isempty(idph) && ~isempty(idmh))
    error('CHEBFUN:DISKFUN:combine:parity',['Inputs must have oposite'...
        'parity. Consider using plus']);
end

pivots = [g.pivotValues(idpg);h.pivotValues(idph);...
          g.pivotValues(idmg);h.pivotValues(idmh)];
cols = [g.cols(:,idpg) h.cols(:,idph) g.cols(:,idmg) h.cols(:,idmh)];
rows = [g.rows(:,idpg) h.rows(:,idph) g.rows(:,idmg) h.rows(:,idmh)];
locations = [g.pivotLocations;h.pivotLocations];

numPlus = length(idpg)+length(idph);
idxPlus = 1:numPlus;
numMinus = length(idmg)+length(idmh);
idxMinus = (numPlus+1):(numPlus+numMinus);

% 
% TODO: The code below sorts columns and rows according to the size of the
% pivots.  This may mess up the assumed location of the non-zero pole term
% for subsequent operations.  Therefore we are commenting it out.  We may
% decide later that sorting is a good idea, in which case we will have to
% figure out how to deal with the non-zero pole.
% Sort the results 
%[ignore,perm] = sort( abs(pivots), 1, 'descend' );
%pivots = pivots(perm);
%cols = cols(:,perm);
%rows = rows(:,perm);

% Figure out where the plus/minus terms went
%idx = 1:length(pivots);
%idx = idx(perm);
%plusFlag = idx <= numel([idpg;idph]);
%idxPlus = find( plusFlag );
%idxMinus = find( ~plusFlag );

% Assemble the results into g.
g.cols = cols;
g.rows = rows;
g.pivotValues = pivots;
g.pivotLocations = locations;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;
g.nonZeroPoles = g.nonZeroPoles || h.nonZeroPoles;

end