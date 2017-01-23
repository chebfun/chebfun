function out = isPeriodicTech(f)
%ISPERIODICTECH   Test if a CHEBFUN2V object is built upon a TECH of 
%periodic functions.
%   out = ISPERIODICTECH(F) returns logical true if the columns and rows of
%   F are made of a TECH of periodic functions and false otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for j = 1:f.nComponents
    out(j) = isPeriodicTech(f.components{j});
end
out = all(out);

end