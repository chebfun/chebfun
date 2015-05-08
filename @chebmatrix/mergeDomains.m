function d = mergeDomains(varargin)
%MERGEDOMAINS  Merge domains of CHEBMATRIX blocks.
%   D = MERGEDOMAINS(D1,D2,...) merges domains (union of breakpoints, while
%   checking endpoints). Each input is inspected.
%       If numeric scalar, it's ignored.
%       If a numeric vector, it's taken as a domain.
%       If a chebfun or linBlock, its domain is extracted.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Start by getting all the domains from VARARGIN
doms = cell(size(varargin));
for argCounter = 1:nargin
    item = varargin{argCounter};
    if ( isnumeric(item) && (length(item) > 1) )
        doms{argCounter} = item;
    elseif ( isa(item, 'chebfun') || isa(item, 'linBlock') )
        doms{argCounter} = item.domain;
    end
end

% Throw away empties
doms(cellfun(@isempty, doms)) = [];

% If we're left with one domain, return it!
if ( length(doms) == 1)
    d = doms{1};
    return
end

% Often in the CHEBMATRIX case, we're dealing with objects that live on the same
% domains. The DOMAIN.MERGE() method is pretty expensive, so we check whether we
% can get away with a cheaper version of merging the domains -- namely, if
% they're all of the same size, we check their diffs, and if they agree, we're
% happy!

% Start off by assuming that the domains don't match:
domainMatch = false;
% 
% Get the lengths of all the domains:
domLengths = cellfun(@length, doms);
% Do the length of the domains agree?
lengthMatch = ~any(domLengths - domLengths(1));
if ( lengthMatch )
    % Compare the other entries of DOMS with the first entry.
    % The row vector with the first domain:
    dom1 = doms{1};
    % A matrix with the rest of the domains:
    domMat = cell2mat(doms(2:end)');
    
    % Subtract the first row from the matrix, and check whether we get any
    % non-zero entries:
    domainMatch = ~any(any(bsxfun(@minus, domMat, dom1)));
    
end

% Did the domains match?
if ( domainMatch )
    % Return the first entry of DOMS:
    d = dom1;
else
    % We need to merge all the domains:
    d = domain.merge(doms{:});
end

end