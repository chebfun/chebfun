function varargout = mldivide(N, rhs, pref, varargin)
%\    MLDIVIDE   Solve CHEBOP BVP or IVP system.
%
% MLDIVIDE is a convenient wrapper for CHEBOP/SOLVEBVP and CHEBOP/SOLVEIVP. It
% calls the appropriate method, depending on whether the problem being solved is
% a boundary-value problem, or an initial/final-value problem. See the
% documentation of the corresponding methods for further details.
%
% See also CHEBOP/SOLVEBVP, CHEBOP/SOLVEIVP.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


% If we did not get preferences passed in, create a cheboppref:
if ( nargin < 3 )
    pref = cheboppref();
end

% Are we dealing with an initial or a final value problem. In that case, either
% we have  
%   - N.LBC is nonempty, but both N.RBC and N.BC are empty. Here, we are dealing
%     with an initial value problem, where all conditions are imposed via N.LBC.
%   - N.RBC is nonempty, but both N.LBC and N.BC are empty. Here, we are dealing
%     with a final value problem, where all conditions are imposed via N.RBC.
isIVP =  ( ~isempty(N.lbc) && isempty(N.rbc) && isempty(N.bc) ) || ...
     ( isempty(N.lbc) && ~isempty(N.rbc) && isempty(N.bc) );
 
% Look at what solver we're using for IVPs:
ivpSolver = func2str(pref.ivpSolver);

% If we are solving an IVP, and we've specified one of the MATLAB solvers
% (ode113, ode15s, ode45) as our choice, call the solveIVP method. Otherwise,
% solve problems globally:
if ( isIVP && ~isempty(strfind(ivpSolver, 'chebfun.ode')) )
    [varargout{1:nargout}] = solveIVP(N, rhs, pref, varargin{:});
else
    % We have conditions in other fields, or we want to solve IVPs globablly,
    % call CHEBOP/SOLVEBVP:
    [varargout{1:nargout}] = solvebvp(N, rhs, pref, varargin{:});
end

end