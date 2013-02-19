function f = fliplr(f)
%FLIPLR  Flip columns of a vectorised FUNCHEB1 object.
%   FLIPLR(F) flips the rows and column of a FUNCHEB1. If F has only one
%   column, then this functon will have no effect.

f.values = fliplr(f.values);
f.coeffs = fliplr(f.coeffs);

end