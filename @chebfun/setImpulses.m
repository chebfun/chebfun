function f = setImpulses(f, j, k, vals)
%SETIMPUSES
%
%   This is needed for Impulse asignment for QUASIMATRICES.
%
%   It may be possible to do something with subsasgn..

% TODO: Document this
% TODO: Deal with array input for quasimatrices
% TODO: Remove this once access to impulses is restricts (when we have deltfun).

if ( numel(f) == 1 )
    % CHEBFUN / array-valued case:
    f.impulses(j,k) = vals;
else
    % QUASIMATRIX case:
    f(j).impulses(k) = vals;
end

end