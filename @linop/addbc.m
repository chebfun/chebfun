function L = addbc(L, varargin)
% Synonym for ADDCONSTRAINT.
%
% See also: linop.addConstraint

L = L.addConstraint(varargin{:});

end