function f = fliplr(f)
%FLIPLR    Flip columns of a vectorised FUN object.
%   FLIPLR(F) flips the columns of a vectorised FUN F in the left/right
%   direction. If F has only one column, then this functon has no effect.

% Flip the ONEFUN:
f.onefun = fliplr(f.onefun);

end