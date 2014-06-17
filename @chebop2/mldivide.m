function varargout = mldivide(varargin)
%\   MLDIVIDE   Solve CHEBOP2 partial differential equations.
%   MLDIVIDE is a convenient wrapper for CHEBOP2/SOLVEPDE, but is limited in
%   that it is adaptive. See CHEBOP2/SOLVEPDE documentation for further 
%   details.
%
% For further details about the PDE solver, see: 
% A. Townsend and S. Olver, The automatic solution of partial differential
% equations using a global spectral method, in preparation, 2014.
% 
% Warning: This PDE solver is an experimental new feature. It has not been
% publicly advertised. 
%
% See also CHEBOP2/SOLVEPDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call CHEBOP2/SOLVEPDE:
[varargout{1:nargout}] = solvepde(varargin{:});

end