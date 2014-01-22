function B = vertcat(varargin)
%VERTCAT   Horizontally concatenate linops.

% Cast the result as a linop if there are any linop arguments.
B = linop( vertcat@chebmatrix(varargin{:}) );

end