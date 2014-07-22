function newDom = merge(varargin)
%MERGE    Merge breakpoints (with a tolerance).
%   MERGE(DOM1, DOM2, ..., DOMN, TOL) merges the breakpoints of the domains
%   DOM1, ..., DOMN to a tolerance TOL, where DOM1, ..., DOMN are sorted vectors
%   of real numbers. If the domains are not compatible, i.e.. the first and
%   final entry of each DOM differ by more than TOL, then an error is returned.
%
%   MERGE(DOM1, DOM2, ..., DOMN) uses a tolerance of the largest magnitude entry
%   in DOM1, ..., DOMN scaled by 10*eps.
%
%   MERGE() always returns an array of doubles, even if the inputs contain
%   DOMAIN objects.
%
% See also WHICHDOMAIN, TWEAKDOMAIN, DOMAINCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Relabel the input variable:
doms = varargin;
% Make sure we have doubles, not domains:
doms = cellfun(@(dom) double(dom), doms, 'UniformOutput', false);

% Ignore empties:
for k = numel(doms):-1:1
    if ( isempty(doms{k}) )
        doms(k) = []; % Discard if this domain is now empty.
    end
end

% Initialise newDom:
newDom = [];

% Choose a tolerance:
if ( isscalar(varargin{end}) )
    % One is given:
    tol = varargin{end};
    doms(end) = [];
else
    % 100*eps*hscale:
    hscales = cellfun(@(d) norm(d, inf), doms);
    hscales(isinf(hscales)) = 1;
    tol = 100*eps*max(hscales);
end
% Check to see if the domains are compatible:
ends = cellfun(@(dom) dom([1 end]), doms, 'UniformOutput', false);
diffEnds = cell2mat(ends.') - repmat(ends{1}, numel(doms), 1);
if ( any(diffEnds(:) > tol) )
    error('CHEBFUN:DOMAIN:merge:incompat', 'Incompatible domains.');
end

j = 1;
while ( ~isempty(doms) )
    % Find the minimum remaining domain entry:
    newDom(j) = min(cellfun(@min, doms)); %#ok<AGROW>
    for k = numel(doms):-1:1
        % Remove -inf
        doms{k}(isinf(doms{k}) & doms{k} < 0) = [];
        % Remove subsequent domains less than tol away:
        doms{k}(doms{k} < newDom(j) + tol) = [];
        if ( isempty(doms{k}) )
            doms(k) = []; % Discard if this domain is now empty.
        end
        if ( isinf(newDom(j)) && newDom(j) > 0 )
            % Stop if we get to +inf
            doms = [];
            break
        end
    end
    j = j + 1;
end 

end
