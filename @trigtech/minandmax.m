function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum on [-1,1].
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the TRIGTECH F on [-1,1]. If F is a
%   array-valued TRIGTECH, VALS is a 2-by-N matrix, where N is the number of
%   columns of F. VALS(1, K) is the global minimum of the Kth column of F on
%   [-1, 1], and VALS(2, K) is the global maximum of the same.
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

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Simply create a CHEBTECH of the TRIGTECH and call minandmax on it.
% We know f will be smooth since it is periodic, thus CHEBTECH should be
% able to beautifully compute the roots.

% An arbitrary decision was made to use CHEBTECH1 in this computation.
g = chebtech1(@(x) f.feval(x));
[vals,pos] = minandmax(g);

end
