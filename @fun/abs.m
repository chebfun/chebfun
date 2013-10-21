function f = abs(f, varargin)
%ABS   Absolute value of a FUN object.
%   ABS(F) returns the absolute value of F, where F is a FUN object with no
%   roots in F.domain. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning. F may be complex.

% Take the absolute value of the ONEFUN:
f.onefun = abs(f.onefun, varargin{:});

end
