function f = fliplr(f)
%FLIPLR   Flip columns of a vectorised FUN object.
%   FLIPLR(F) flips the columns of a vectorised FUN in the left/right
%   direction. If F has only one column, then this functon has no effect.

f.onefun = fliplr(f.onefun);

end