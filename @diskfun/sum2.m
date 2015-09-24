function v = sum2( f ) 
% Definite integration of a diskfun. 

% Split f into its plus/minus terms.  The minus terms have integral zero 
% on the disk since the rows in this case are anti-periodic with period
% pi.  Thus, we only need to integrate the plus terms.
[cols,d,rows] = cdr(f);

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

%
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