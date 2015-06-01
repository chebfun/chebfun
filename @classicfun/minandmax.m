function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum of a CLASSICFUN on [a,b].
%   VALS = MINANDMAX(F) returns a 2-vector, VALS = [MIN(F); MAX(F)], with the
%   global minimum and maximum of the CLASSICFUN F on [a,b].  If F is an array-valued
%   CLASSICFUN then VALS is a 2-by-M matrix, where M is the number of columns of F.
%   VALS(1, K) is the global minimum of the Kth column of F on [a,b], and
%   VALS(2, K) is the global maximum of the same.
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and the maximum of F occur.
%
%   If F is complex-valued the absolute values are taken to determine extrema
%   but the output values correspond to those of the original function. That is,
%   VALS = FEVAL(F, POS) where [~, POS] = MINANDMAX(ABS(F)).
%
% See also MIN, MAX.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call minandmax for the ONEFUN of f:
[vals, onefunPos] = minandmax(f.onefun);

% Map the position results obtained on [-1,1] to [a,b]:
pos = f.mapping.For(onefunPos);

end
