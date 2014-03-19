function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum on [-pi,pi].
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the FOURIERTECH F on [-1,1].  If F is a
%   array-valued FOURIERTECH, VALS is a 2-by-N matrix, where N is the number of
%   columns of F.  VALS(1, K) is the global minimum of the Kth column of F on
%   [-pi, pi], and VALS(2, K) is the global maximum of the same.
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and maximum of F occur.
%
%   If F is complex-valued the absolute values are taken to determine extrema
%   but the resulting values correspond to those of the original function. That
%   is, VALS = FEVAL(F, POS) where [~, POS] = MINANDMAX(ABS(F)). (In fact,
%   MINANDMAX actually computes [~, POS] = MINANDMAX(ABS(F).^2), to avoid
%   introducing singularities to the function).
%
% See also MIN, MAX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Simply create a chebfun of the FOURIERTECH and call minandmax on it.
g = chebfun(@(x) f.feval(x),[-pi,pi]);
[vals,pos] = minandmax(g);

end
