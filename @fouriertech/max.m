function [maxVal, maxPos] = max(f)
%MAX   Global maximum of a FOURIERTECH on [-pi,pi].
%   MAXVAL = MAX(F) returns the global maximum of the FOURIERTECH on [-pi,pi].  If
%   F is an array-valued FOURIERTECH, MAXVAL is a row vector whose Kth entry is the
%   global maximum of the Kth column of F.
%
%   [MAXVAL, MAXPOS] = MAX(F) returns also a value such that MAXVAL = F(MAXPOS).
%
%   If F is complex-valued then absolute values are taken to determine maxima
%   but the resulting value corresponds to that of the original function. That
%   is, MAXVAL = feval(F, MAXPOS) where [~, MAXPOS] = MAX(abs(F));
%
% See also MIN, MINANDMAX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Simply create a chebfun of the FOURIERTECH and call max on it.
g = chebfun(@(x) f.feval(x),[-pi,pi]);
[maxVal, maxPos] = max(g);

end
