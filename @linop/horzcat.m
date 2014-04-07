function B = horzcat(varargin)
%HORZCAT   Horizontally concatenate linops.

% Cast the result as a linop if there are any linop arguments.
B = linop( horzcat@chebmatrix(varargin{:}) );

end