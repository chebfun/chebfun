function varargout = bvp5c(fun1, fun2, y0, varargin)
%BVP5C   Solve boundary value problems for ODEs by collocation with CHEBFUN.
%   Y = BVP5C(ODEFUN, BCFUN, Y0) applies the standard BVP5C method to solve a
%   boundary-value problem. ODEFUN and BCFUN are as in BVP5C. The Y0 argument is
%   a CHEBFUN that represents the initial guess to the solution Y. Its domain
%   defines the domain of the problem, and the length of the CHEBFUN Y0 is used
%   to set the number of points in an initial equispaced mesh. Note that it is
%   not necessary to call BVPINIT.
%
%   [Y, P] = BVP5C(ODEFUN, BCFUN, Y0, PARAM, OPTS) allows you to specify an
%   initial guess for any additional parameters to be found for the solution,
%   and an options vector to guide the solution. See the built in BVP5C and
%   BVPSET for details. You may specify either extra argument, or both. An
%   additional output is used to return the parameter values found.
%
%   It is possible to take a crude continuation approach by solving for a simple
%   variation of the problem, then using the resulting CHEBFUN as the initial
%   guess for a more difficult version.
%
%   Note that CHEBFUN/BVP5C() uses a default RELTOL of 1e-6.
%
% Example (using built-in BVP demo functions):
%   y0 = chebfun([0, 0], [0, 4]);
%   y = bvp5c(@twoode, @twobc, y0);
%   plot(y)
%
% See also BVPINIT, BVPSET, BVP4C, ODE113.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse the inputs.
params = {};
opts = {};
for k = 1:nargin-3
    t = varargin{k};
    if ( isstruct(t) )
        opts{1} = t;
    elseif ( isnumeric(t) )
        params{1} = t;
    end
end

if ( ~isfinite(y0) )
    error('CHEBFUN:CHEBFUN:bvp5c:inf',...
      'BVP5C() does not currently support functions which diverge to infinity');
end

% Determine the initial BVP grid from the length of the CHEBFUN:
n = max(9, length(y0));
dom = domain(y0);
x = linspace(dom(1), dom(end), n);

% Call BVPINT().
f = @(x) feval(y0, x);
init = bvpinit(x, f, params{:});

% Call BVP solver and convert to CHEBFUN.
sol = bvp5c(fun1, fun2, init, opts{:});
varargout{1} = chebfun.odesol(sol, dom, opts{:});

% Look for parameter output.
if ( ~isempty(params) )
    varargout{2} = sol.parameters; 
end

end
