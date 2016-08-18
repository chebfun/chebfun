function out = isPeriodicTech(f)
%ISPERIODICTECH   Test if a CHEBFUN2 object is built upon a TECH of periodic
%functions.
%   out = ISPERIODICTECH(F) returns logical true if the columns and rows of F
%   are made of a TECH of periodic functions and false otherwise.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isPeriodicTech(f.cols) .* isPeriodicTech(f.rows);

% TO DO: would the following be safer?
% out = 1;
% for jj = 1:rank(f)
%     out = out .* isPeriodicTech(f.cols(:,jj));
%     out = out .* isPeriodicTech(f.rows(:,jj));
% end

end
