function B = cat(varargin)

% Cast the result as a linop if there are any linop arguments.

B = linop( cat@chebmatrix(varargin{:}) );

end