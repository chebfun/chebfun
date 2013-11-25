function A = discretize(L, varargin)

% Construct a single matrix based on DIM-dimensional blocks.

if isa(varargin{1},'blockDiscretization')
    dsc = varargin{1};
else
    dsc = L.discretizationType( varargin{:} );
end

% The chebmatrix domain overrides that in the given discretization object, if
% any.
dsc.domain = L.domain;

A = cellfun(@(x) matrix(dsc,x),A.blocks,'uniform',false);

end
