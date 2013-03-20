function f = fliplr(f)
%FLIPLR   Flip columns of a vectorised CHEBTECH object.
%   FLIPLR(F) flips the columns of a vectorised CHEBTECH in the left/right
%   direction. If F has only one column, then this functon has no effect.

f.values = fliplr(f.values);
f.coeffs = fliplr(f.coeffs);

end
