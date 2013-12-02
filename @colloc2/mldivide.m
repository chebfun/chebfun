function [v,disc] = mldivide(disc,A,b)
s = 1./ max(1, max(abs(A),[],2) );
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
A = bsxfun(@times,s,A);
v = A \ (s.*b);
end
