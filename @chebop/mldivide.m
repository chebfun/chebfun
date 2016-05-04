function varargout = mldivide(N, rhs, varargin)
%\    MLDIVIDE   Solve CHEBOP BVP or IVP system.
%
% MLDIVIDE is a convenient wrapper for CHEBOP/SOLVEBVP and CHEBOP/SOLVEIVP. It
% calls the appropriate method, depending on whether the problem being solved is
% a boundary-value problem, or an initial/final-value problem. Problems,
% specified by a CHEBOP N are determined to be a boundary-value or
% initial/final-value problems as follows:
%   * If N.LBC is non-empty, but both N.RBC and N.BC are empty, it's considered
%     to be an initial-value problem. In this case, CHEBOP/SOLVEIVP is called.
%   * If N.RBC is non-empty, but both N.LBC and N.BC are empty, it's considered
%     to be a final-value problem. In this case, CHEBOP/SOLVEIVP is called.
%   * Otherwise, the problem is considered to be a boundary-value problem. In
%     this case, CHEBOP/SOLVEBVP is called.
%
% See the documentation of the corresponding methods for further details.
%
% See also CHEBOP/SOLVEBVP, CHEBOP/SOLVEIVP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Are we dealing with an initial or a final value problem. In that case, either
% we have  
%   - N.LBC is nonempty, but both N.RBC and N.BC are empty. Here, we are dealing
%     with an initial value problem, where all conditions are imposed via N.LBC.
%   - N.RBC is nonempty, but both N.LBC and N.BC are empty. Here, we are dealing
%     with a final value problem, where all conditions are imposed via N.RBC.
isIVP =  ( ~isempty(N.lbc) && isempty(N.rbc) && isempty(N.bc) ) || ...
     ( isempty(N.lbc) && ~isempty(N.rbc) && isempty(N.bc) );

% If we did not get preferences passed in, we need to create a cheboppref so
% that we know wheter we want to solve the IVP via time-stepping or globally:
if ( nargin < 3 )
    pref = cheboppref();
else
    pref = varargin{1};
end

% Look at what solver we're using for IVPs:
ivpSolver = func2str(pref.ivpSolver);

% If we are solving an IVP, and we've specified one of the MATLAB solvers
% (ode113, ode15s, ode45) as our choice, call the solveivp method. Otherwise,
% solve problems globally:
if ( isIVP && ~isempty(strfind(ivpSolver, 'chebfun.ode')) )
    % Pass PREF here, as we'd have to start by constructing a CHEBOPPREF anyway
    % that start of solveivp otherwise.
    [varargout{1:nargout}] = solveivp(N, rhs, pref, varargin{:});
else
    % We have conditions in other fields, or we want to solve IVPs globablly,
    % call CHEBOP/SOLVEBVP. However, here, we don't want to pass PREF if it
    % wasn't included in VARARGIN, as that'd mean we couldn't to an automatic
    % switch to Fourier methods in the periodic case (i.e. change the
    % discretization to TRIGTECH).
    [varargout{1:nargout}] = solvebvp(N, rhs, varargin{:});
end

end
