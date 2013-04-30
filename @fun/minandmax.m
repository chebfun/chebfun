function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum of a FUN on [a,b].
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the FUN F on [a,b].  If F is an array-valued
%   FUN, VALS is a 2-by-N matrix, where N is the number of columns of F.
%   VALS(1, K) is the global minimum of the Kth column of F on [a,b], and
%   VALS(2, K) is the global maximum of the same.
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and the maximum of F occur.
%
%   If F is complex-valued the absolute values are taken to determine extrema
%   but the output values correspond to those of the original function. That is,
%   VALS = FEVAL(F, POS) where [~, POS] = MINANDMAX(ABS(F)). (In fact, MINANDMAX
%   actually computes [~, POS] = MINANDMAX(ABS(F).^2), to avoid introducing
%   singularlities to the function).
%
% See also MIN, MAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call minandmax for the onefun of f:
[vals, onefunPos] = minandmax(f.onefun);

% Map the position results obtained on [-1,1] to [a,b]:
pos = f.mapping.for(onefunPos);

end
