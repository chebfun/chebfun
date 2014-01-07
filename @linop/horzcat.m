function B = horzcat(varargin)

% Cast the result as a linop if there are any linop arguments.

B = linop( horzcat@chebmatrix(varargin{:}) );

end