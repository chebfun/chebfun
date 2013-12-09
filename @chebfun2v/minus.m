function f = minus(f,g)
% - MINUS. Minus of two chebfun2v.  
%
% f - g substracts the chebfun2v f from g componentwise. 
% minus(f,g) is called for the syntax f - g. 

% Copyright 2012 by The University of Oxford and The Chebfun2 Developers. 

f = plus(f,uminus(g));

end