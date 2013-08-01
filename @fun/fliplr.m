function f = fliplr(f)
%FLIPLR    Flip columns of an array-valued FUN object.
%   FLIPLR(F) flips the columns of an array-valued FUN F in the left/right
%   direction. If F has only one column, then this function has no effect.

% Flip the ONEFUN:
f.onefun = fliplr(f.onefun);

end
