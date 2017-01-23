function out = isPeriodicTech(f)
%ISPERIODICTECH   Test if a CHEBFUN2 object is built upon a TECH of periodic
%functions.
%   out = ISPERIODICTECH(F) returns logical true if the columns and rows of F
%   are made of a TECH of periodic functions and false otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isPeriodicTech(f.cols) & isPeriodicTech(f.rows);

end