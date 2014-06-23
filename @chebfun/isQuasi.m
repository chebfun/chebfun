function out = isQuasi(f)
%ISQUASI   Returns true for a quasimatrix input and flase otherwise.
%   ISQUASI(F) returns true for a quasimatrix input and false otherwise. This is
%   essentially a wrapper for NUMEL(F) > 1, but is convenient for explaining how
%   quasimatrices work.
%
% See also QUASIMATRIX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = numel(f) > 1;

end