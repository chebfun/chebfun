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

s = cellfun(@(b) min(size(b)), A.blocks);

% If A contains inf x inf blocks, call SPY().
if ( ~all(isfinite(s(:))) )
    spy(A);
    
% If A contains only CHEBFUN or DOUBLE, convert it to a QUASIMATRIX, and
% cal CHEBFUN/PLOT().
else
    A = chebfun(A);
    plot(A);
end

end