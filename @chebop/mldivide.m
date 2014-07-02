function varargout = mldivide(N, varargin)
%\    MLDIVIDE   Solve CHEBOP BVP system.
%   MLDIVIDE is a convenient wrapper for CHEBOP/SOLVEBVP, but is limited in that
%   it only supports a single output. See CHEBOP/SOLVEBVP documentation for
%   further details.
%
% See also CHEBOP/SOLVEBVP.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(N.rbc) && isempty(N.bc) )
    % Both N.RBC and N.BC are empty, we must be dealing with an IVP where all
    % conditions are imposed via N.LBC
    [varargout{1:nargout}] = solveivp(N, varargin{:});
else
    % We have conditions in other fields, call CHEBOP/SOLVEBVP:
    [varargout{1:nargout}] = solvebvp(N, varargin{:});
end

end