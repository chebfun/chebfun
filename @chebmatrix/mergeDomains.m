function d = mergeDomains(varargin)
%MERGEDOMAINS  Merge domains of CHEBMATRIX blocks.
%   D = MERGEDOMAINS(D1,D2,...) merges domains (union of breakpoints, while
%   checking endpoints). Each input is inspected.
%       If numeric scalar, it's ignored.
%       If a numeric vector, it's taken as a domain.
%       If a chebfun or linBlock, its domain is extracted.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

d = [];
for i = 1:nargin
    item = varargin{i};
    if ( isnumeric(item) && (length(item) > 1) )
        d = domain.merge(d, item);
    elseif ( isa(item, 'chebfun') || isa(item, 'linBlock') )
        d = domain.merge(d, item.domain);
    end
end

end