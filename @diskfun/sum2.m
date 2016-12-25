function v = sum2( f ) 
%SUM2   Double integral of a DISKFUN over its domain.
%   I = SUM2(F) returns the double definite integral of a DISKFUN.


% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[cols,d,rows] = cdr(f);

% Split f into its plus/minus terms.  The minus terms have integral zero 
% on the disk since the rows in this case are anti-periodic with period
% pi. Thus, we only need to integrate the plus terms.

% If there are no plus terms then the integral is zero
if isempty(f.idxPlus)
    v = 0;
    return
end

cols = cols(:,f.idxPlus);
rows = rows(:,f.idxPlus);
d = diag(d(f.idxPlus,f.idxPlus)).';

% Integrate the rows over their domain.
intRows = sum(rows);

% Create a chebfun of the measure. 
measure = chebfun(@(r) r,[-1,1]);

% Multiply the columns by the measure
cols = cols.*(measure*ones(1,size(cols,2)));

% Integrate each column over the non-doubled up r coordinate.
intCols = sum(cols,[0 1]);

% Put the integrals together to get the final result.
v = sum(d.*intRows.*intCols);

% TODO: Add support for different domains.

end