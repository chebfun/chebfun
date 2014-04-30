function plot(A)
%PLOT   Plot for CHEBMATRIX objects.
%   PLOT(A) plots the CHEBMATRIX object A.
%
%   If A contains only CHEBFUN and DOUBLE objects, A is converted to a
%   QUASIMATRIX, and CHEBFUN/PLOT is called.
%
%   If not, CHEBMATRIX/SPY is called.
%
%   See also CHEBMATRIX, CHEMATRIX/SPY, CHEBFUN/PLOT.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

sz = size(A, 1)*size(A, 2);
temp = 1;

% Check if A contains only CHEBFUN and DOUBLE objects.
for j = 1:sz
   if  ( isa(A.blocks{j}, 'chebfun') | isa(A.blocks{j}, 'double') )
   else
      temp = 0;    
   end
end

% If so, convert A to a QUASIMATRIX.
if temp == 1
   F = chebfun(A);
   plot(F);

% If not, cal CHEBMATRIX/SPY.
else
   spy(A);
end

end