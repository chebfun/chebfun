function d = mergeDomains(A,varargin)

% D = MERGEDOMAINS(A,D1,D2,...), where D1,... are numeric vectors, performs the
% result of a chebfun.mergeDomains on the numeric domains with that of the
% chebmatrix A.
%
% D = MERGEDOMAINS(A1,A2,...) does the same using the domains of the
% chebmatrices A1,A2,....

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin == 1 )
    d = A.domain;
elseif isnumeric(varargin{1})
    d = chebfun.mergeDomains(A.domain,varargin{:});
else
    c = cellfun(@(x) x.domain,varargin,'uniform',false);
    d = chebfun.mergeDomains(c{:});
end

end
