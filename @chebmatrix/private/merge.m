function d = merge(varargin)
% D = MERGEDOMAINS(D1,D2,...) merges domains (union of breakpoints,
% while checking endpoints). Each input is inspected.
%   If numeric scalar, it's ignored.
%   If a numeric vector, it's taken as a domain.
%   If a chebfun or linBlock, its domain is extracted.

d = [];
for i = 1:nargin
    item = varargin{i};
    if ( isnumeric(item) && length(item) > 1 )
        d = chebfun.mergeDomains(d,item);
    elseif ( isa(item,'chebfun') || isa(item,'linBlock') )
        d = chebfun.mergeDomains(d,item.domain);
    end
end
end