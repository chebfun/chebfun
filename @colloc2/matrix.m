function varargout = matrix(dsc,dimension,domain)
if nargin > 1
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    dsc.dimension = dimension;
    if nargin > 2
        dsc.domain = domain;
    end
end
L = dsc.source;
if isa(L,'chebmatrix')
    A = cellfun(@(x) blockMatrix(dsc,x),L.blocks,'uniform',false);
    out{1} = A;
    if isa(L,'linop')
        [out{1:3}] = useConstraints(dsc,A);
    end
    m = max(1,nargout);
    varargout(1:m) = out(1:m);
else
    [varargout{1:nargout}] = blockMatrix(dsc);
end
end
