function [normA, normLoc] = norm(A, n)
%NORM   Norm of a CHEBMATRIX object.
%   NORM(A) computes the norm of the CHEBMATRIX object A.
%
%   If A contains only CHEBFUN and DOUBLE objects, A is converted to a
%   QUASIMATRIX, and CHEBFUN/NORM is called.
%
%   If not, CHEBMATRIX/? is called. [TODO]
%
%   See also CHEBMATRIX, CHEBFUN/NORM.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Empty CHEBMATRIX has norm 0:
if ( isempty(A) )
    normA = 0;
    return
end

if ( nargin == 1 )
    n = 'fro'; 	% Frobenius norm is the default.
end

sz = size(A, 1)*size(A, 2);
temp = 1;

% Check if A contains only CHEBFUN and DOUBLE objects.
for j = 1:sz
   if  ( isa(A.blocks{j}, 'chebfun') | isa(A.blocks{j}, 'double') )
   else
      temp = 0;    
   end
end

% If so, convert A to a QUASIMATRIX, and call CHEBFUN/NORM.
if temp == 1
   F = chebfun(A);
   [normA, normLoc] = norm(F, n);

% If not, ?. [TODO]
else
   normA = 0;
   normLoc = 0;
end

end