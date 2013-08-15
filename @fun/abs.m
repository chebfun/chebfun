function f = abs(f)
%ABS   Absolute value of a FUN object.
%   % ABS(F) returns the absolute values of F, where F is a FUN object with no
%   roots in F.domain. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning.

% Take the absolute value of the ONEFUN:
f.onefun = abs(f.onefun);

end