function B = vertcat(varargin)
%VERTCAT   Horizontally concatenate LINOP objects.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Cast the result as a linop if there are any linop arguments.
B = linop( vertcat@chebmatrix(varargin{:}) );

end
