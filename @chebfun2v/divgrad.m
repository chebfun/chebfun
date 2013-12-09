function f = divgrad(f)
%DIVGRAD Laplacian of a chebfun2v.
%
% F = DIVGRAD(F) returns the Laplacian of a chebfun2v i.e.
% 
% divgrad(f) = f_xx + f_yy 
%
% Also see chebfun2v/lap

% Copyright 2012 by The University of Oxford and The Chebfun2 Developers. 

f = diff(f.xcheb,[2,0]) + diff(f.ycheb,[0,2]); % laplacian. 
 
end