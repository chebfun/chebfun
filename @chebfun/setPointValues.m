function f = setPointValues(f, j, k, vals)
%SETPOINTVALUES
%
%   This is needed for Impulse asignment for QUASIMATRICES.
%
%   It may be possible to do something with subsasgn..

% TODO: Document this
% TODO: Deal with array input for quasimatrices
% TODO: Remove this once access to impulses is restricts (when we have deltfun).

if ( numel(f) == 1 )
    % CHEBFUN / array-valued case:
    f.pointValues(j,k) = vals;
else
    % QUASIMATRIX case:
    f(j).pointValues(k) = vals;
end

end