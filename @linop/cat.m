function B = cat(varargin)
%CAT    Concatenate LINOP objects.

% Cast the result as a linop if there are any linop arguments.
B = linop( cat@chebmatrix(varargin{:}) );

end