function f = fliplr(f)
%FLIPLR  Flip columns of a vectorised FUNCHEB2 object.
%   FLIPLR(F) flips the columns of a vectorised FUNCHEB2 in the left/right
%   direction. If F has only one column, then this functon will have no effect.

f.values = fliplr(f.values);
f.coeffs = fliplr(f.coeffs);

end
