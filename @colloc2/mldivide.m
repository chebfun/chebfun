function [v,disc] = mldivide(disc,A,b)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Can we reuse previously computed LU factors?


if ( ~isFactored(disc) )
    s = 1./ max(1, max(abs(A),[],2) );
    A = bsxfun(@times,s,A);
    [L,U] = lu(A);
    disc.mldivideData = {L,U,s};
else
    [L,U,s] = deal(disc.mldivideData{:});
end

v = U \ ( L\ (s.*b) );

end
