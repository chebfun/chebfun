classdef chebfun
%CHEBFUN   CHEBFUN class for representing functions on [a,b].
%
%   Class for approximating functions defined on finite, semi-infinite, or
%   doubly-infinite intervals [a,b]. Functions may be smooth, piecewise smooth,
%   weakly singular, or blow up on the interval.
%
% CHEBFUN(F) constructs a CHEBFUN object representing the function F on the
% interval [-1,1]. F may be a string, e.g., 'sin(x)', a function handle, e.g.,
% @(x) x.^2 + 2*x + 1, or a vector of numbers. In the first two instances, F
% should be "vectorized" in the sense that it may be evaluated at a column
% vector of points x(:) in [-1,1] and return an output of size NxM where N =
% length(x(:)). If this is not possible then the flag CHEBFUN(F, 'vectorize')
% should be passed. CHEBFUN(F, 'vectorcheck', 'off') disables the automatic
% checking for vector input. CHEBFUN() returns an empty CHEBFUN object.
%
% CHEBFUN(F, [A, B]) specifies an interval [A,B] on which the CHEBFUN is
% defined, where A and/or B may be infinite. CHEBFUN(F, ENDS), where ENDS is a
% 1x(K+1) vector of unique ascending values, specifies a piecewise smooth
% CHEBFUN defined on the interval [ENDS(1), ENDS(K+1)] with additional interior
% breaks at ENDS(2), ..., ENDS(K). Specifying these breaks can be particularly
% useful if F is known to have discontinuities. For example,
%   CHEBFUN(@(x) abs(x), [-1, 0, 1]).
% If a domain is passed to the constructor, it should always be the 2nd input.
%
% CHEBFUN(A) or CHEBFUN(A, 'chebpts2'), where A is an Nx1 matrix, constructs a
% CHEBFUN object which interpolates the data in A on an N-point Chebyshev grid
% of the second kind (see >> help chebpts). CHEBFUN(A, 'chebpts1') and
% CHEBFUN(A, 'equi') are similar, but here the data is assumed to come from a
% 1st-kind Chebyshev or equispaced grid linspace(-1, 1, N), respectively. (In
% the latter case, a smooth interpolant is constructed using an adaptive
% Floater-Hormann scheme [Numer. Math. 107, 315-331 (2007)].). CHEBFUN(F, N) or
% CHEBFUN(F, N, 'chebpts2') is equivalent to CHEBFUN(feval(F, chebpts(N)).
%
% CHEBFUN({F1,...,Fk}, ENDS) constructs a piecewise smooth CHEBFUN which
% represents Fj on the interval [ENDS(j), END(j+1)]. Each entry Fj may be a
% string, function handle, or vector of doubles. For example
%   CHEBFUN({@(x) sin(x), @(x) cos(x)}, [-1, 0, 1])
%
% CHEBFUN(F, PREF) or CHEBFUN(F, [A, B], PREF) constructs a CHEBFUN object from
% F with the options determined by the CHEBPREF object PREF.
% Construction time options may also be passed directly to the constructor in
% the form CHEBFUN(F, [A, B], PROP1, VAL1, PROP2, VAL2, ...). (See
% CHEBPREF for details of the various preference options.). In particular,
% CHEBFUN(F, 'splitting', 'on') allows the constructor to adaptively determine
% breakpoints to better represent piecewise smooth functions F. For example,
%   CHEBFUN(@(x) sign(x - .3), [-1, 1], 'splitting', 'on')
% It is not possible to mix PROP/VAL and PREF inputs in a single constructor
% call.
%
% CHEBFUN(F, ...), where F is an NxM matrix or an array-valued function handle,
% returns an "array-valued" CHEBFUN. For example,
%   CHEBFUN(rand(14, 2))
% or
%   CHEBFUN(@(x) [sin(x), cos(x)])
% Note that each column in an array-valued CHEBFUN object is discretized in the
% same way (i.e., the same breakpoint locations and the same underlying
% representation). Note the difference between 
%   CHEBFUN(@(x) [sin(x), cos(x)], [-1, 0, 1])
% and
%   CHEBFUN({@(x) sin(x), @(x) cos(x)}, [-1, 0, 1]).
% The former constructs an array-valued CHEBFUN with both columns defined on the
% domain [-1, 0, 1]. The latter defines a single column CHEBFUN which represents
% sin(x) in the interval [-1, 0) and cos(x) on the interval (0, 1].
%
% See also CHEBPREF, CHEBPTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBFUN Class Description:
%
% The CHEBFUN class is for representations of piecewise functions on the
% interval [a,b].
%
% The CHEBFUN class is the main user interface. We do not expect users to
% directly invoke any objects below this level.
%
% A CHEBFUN object consists of a collection of FUN objects. There are two main
% tasks for the CHEBFUN constructor: (1) parse the user input, and (2) correctly
% piece together FUN objects to form a global approximation. If the input
% function is globally smooth then the resulting CHEBFUN contains a single FUN
% object. If the input is not smooth, or breakpoints are passed to the
% constructor, CHEBFUN must determine appropriate breakpoints and return a
% piecewise smooth CHEBFUN with multiple FUN objects.
%
% This is a user-level class, and all input arguments should be thoroughly
% sanity checked.
%
% Class diagram: [ADchebfun] <>-- [CHEBFUN] <>-- [<<fun>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Access = public)
        % DOMAIN of definition of a CHEBFUN object. If K = length(F.DOMAIN) is
        % greater than 1 then the CHEBFUN is referred to as a "piecewise".
        % CHEBFUN. The first and last values of this vector define the left and
        % right endpoints of the domain, respectively. The other values give the
        % locations of the interior breakpoints that define the domains of the
        % individual FUN objects comprising the CHEBFUN. The entries in this
        % vector should be strictly increasing.
        domain              % (1x(K+1) double)

        % FUNS is a cell array containing the FUN objects that comprise a
        % piecewise CHEBFUN. The the kth entry in this cell is the FUN defining
        % the representation used by the CHEBFUN object on the open interval
        % (F.DOMAIN(k), F.DOMAIN(k+1)). If M = size(f.funs, 2) is greater than
        % 1, then the CHEBFUN object is referred to as "array valued".
        funs                % (Kx1 cell array of FUN objects)

        % IMPULSES is a three-dimensional array storing information about the
        % values of the CHEBFUN object at the points in DOMAIN. The rows
        % correspond to the breakpoints in the DOMAIN vector, and if M > 1 then
        % the columns correspond to the columns in an array-valued CHEBFUN.
        % Thus, F.IMPULSES(:, :, 1) is a matrix consisting of the values of
        % each column of F at each breakpoint. The third dimension is used for
        % storing information about higher-order delta functions that may be
        % present at breakpoints. (See "help dirac" for more details.)
        impulses = [];      % ((K+1) x M x (D+1) double)

        % ISTRANSPOSED determines whether a (possibly array-valued) CHEBFUN F
        % should be interpreted as a collection of "column" CHEBFUN objects (if
        % F.ISTRANSPOSED == 0, the default), which are considered (infxM)
        % arrays, or "row" CHEBFUN objects (if F.ISTRANSPOSED == 1), which are
        % (Mxinf) arrays. This difference is only behavioral; the other
        % properties described above are _NOT_ stored differently if this flag
        % is set.)
        isTransposed = 0;   % (logical)
    end
    
    methods
        
        function f = chebfun(varargin)
            % The main CHEBFUN constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            % Assign variable to workspace for string input:
            % TODO: Document this.
            if ( all(cellfun(@isvarname, varargin)) )
                f = {};
                for k = numel(varargin):-1:1
                    f{k} = chebfun(str2op(varargin{k}));
                    assignin('caller', varargin{k}, f{k});
                end
                f = horzcat(f{:});
                return
            end
            
            % Parse inputs:
            [op, dom, pref] = parseInputs(varargin{:});
            
            if ( iscell(op) && all(cellfun(@(x) isa(x, 'fun'), op)) )
                % Construct a CHEBFUN from a cell array of FUN objects:
                
                if ( nargin > 1 )
                    error('CHEBFUN:chebfun:nargin', ...
                        'Only one input is allowed when passing an array of FUNs.')
                end
                
                % Assign the cell to the .FUNS property:
                f.funs = op;
                % Collect the domains together:
                dom = cellfun(@(fun) get(fun, 'domain'), f.funs, ...
                    'uniformOutput', false);
                f.domain = unique([dom{:}]);
                % Update values at breakpoints (first row of f.impulses):
                f.impulses = chebfun.getValuesAtBreakpoints(f.funs, f.domain);
                
            else
                % Construct from function_handle, numeric, or string input:
                
                % Call the main constructor:
                [f.funs, f.domain] = chebfun.constructor(op, dom, pref);
                
                % Update values at breakpoints (first row of f.impulses):
                f.impulses = chebfun.getValuesAtBreakpoints(f.funs, f.domain, op);
                
                % Remove unnecessary breaks (but not those that were given):
                [ignored, index] = setdiff(f.domain, dom);
                f = merge(f, index(:).', pref);
                
            end
            
        end
        
    end
    
    % Static methods implemented by CHEBFUN class.
    methods (Static = true)
        
        % Main constructor.
        [funs, ends] = constructor(op, domain, pref);
        
        % Edge detector.
        [edge, vscale] = detectEdge(op, domain, hscale, vscale, derHandle);
        
        % Interpolate data:
        f = interp1(x, y, method, dom);
        
        % Determine values of chebfun at breakpoints.
        vals = getValuesAtBreakpoints(funs, ends, op);
        
        % ODE113 with CHEBFUN output.
        [t, y] = ode113(varargin);
        
        % ODE15S with CHEBFUN output.
        [t, y] = ode15s(varargin);
        
        % ODE45 with CHEBFUN output.
        [t, y] = ode45(varargin);
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin);
        
        % Cubic Hermite interpolation:
        f = pchip(x, y, method);
        
        % Cubic spline interpolant:
        f = spline(x, y, d);
        
    end

    methods (Access = private)
        % Remove zero layers from impulses array.
        f = tidyImpulses(f);

        % Set small breakpoint values to zero.
        f = thresholdBreakpointValues(f);
    end
    
    % Static private methods implemented by CHEBFUN class.
    methods (Static = true, Access = private)
        
        % Convert ODE solutions into CHEBFUN objects:
        [y, t] = odesol(sol, opt);
        
        % Parse the inputs to the CHEBFUN constructor.
        [op, domain, pref] = parseInputs(op, domain, varargin);
        
        % Convert a string input to a function_handle.
        op = str2op(op);
        
        % Vectorise a function handle input.
        op = vec(op);
        
    end
    
    % Methods implemented by CHEBFUN class.
    methods
        % True if any element of a CHEBFUN is a nonzero number, ignoring NaN.
        a = any(f, dim)

        % Absolute value of a CHEBFUN.
        f = abs(f, pref)

        % Round a CHEBFUN towards plus infinity.
        g = ceil(f)

        % Solve boundary value problems for ODEs by collocation.
        [y, t] = bvp4c(fun1, fun2, y0, varargin);

        % Solve boundary value problems for ODEs by collocation.
        [y, t] = bvp5c(fun1, fun2, y0, varargin);

        % Plot information regarding the representation of a CHEBFUN object:
        h = chebpolyplot(f, varargin);

        % Construct complex CHEBFUN from real and imaginary parts.
        C = complex(A, B)

        % Compose CHEBFUN objects with another function.
        h = compose(f, op, g, pref)

        % Compose two CHEBFUN objects (i.e., f(g)).
        h = composeChebfuns(f, g, pref)
        
        % Complex conjugate of a CHEBFUN.
        f = conj(f)
        
        % Complex transpose of a CHEBFUN.
        f = ctranspose(f)

        % Display a CHEBFUN object.
        display(f);
        
        % Compare domains of two CHEBFUN objects.
        pass = domainCheck(f, g);

        % Accuracy estimate of a CHEBFUN object.
        out = epslevel(f);
        
        % Extract columns of an array-valued CHEBFUN object.
        f = extractColumns(f, columnIndex);
        
        % Evaluate a CHEBFUN.
        y = feval(f, x, varargin)
        
        % Round a CHEBFUN towards zero.
        g = fix(f);
        
        % Round a CHEBFUN towards minus infinity.
        g = floor(f);

        % Get properties of a CHEBFUN object.
        out = get(f, prop);
        
        % Horizontal scale of a CHEBFUN object.
        out = hscale(f);

        % Imaginary part of a CHEBFUN.
        f = imag(f)
        
        % True for an empty CHEBFUN.
        out = isempty(f)

        % Test if CHEBFUN objects are equal.
        out = isequal(f, g)

        % Test if a CHEBFUN is bounded.
        out = isfinite(f)
        
        % Test if a CHEBFUN is unbounded.
        out = isinf(f)

        % Test if a CHEBFUN has any NaN values.
        out = isnan(f)
        
        % True for real CHEBFUN.
        out = isreal(f);
        
        % True for zero CHEBFUN objects
        out = iszero(f)
        
        % Length of a CHEBFUN.
        out = length(f);
        
        % Return Legendre coefficients of a CHEBFUN.
        c_leg = legpoly(f, n)
        
        % Plot a CHEBFUN object on a loglog scale:
        h = loglog(f, varargin);
        
        % Plot a CHEBFUN object:
        varargout = plot(f, varargin);
        
        % 3-D plot for CHEBFUN objects.
        varargout = plot3(f, g, h, varargin)
        
        % Subtraction of two CHEBFUN objects.
        f = minus(f, g)
        
        % Signmum of a CHEBFUN.
        f = sign(f, pref)

        % Multiplication of CHEBFUN objects.
        f = mtimes(f, c)

        % Remove unnecessary breakpoints in from a CHEBFUN.
        [f, mergedPts] = merge(f, index, pref)
        
        % Overlap the domain of two CHEBFUN objects.
        [f, g] = overlap(f, g)

        % Obtain data used for plotting a CHEBFUN object:
        data = plotData(f, g, h)
        
        % Power of a CHEBFUN
        f = power(f, b);
        
        % Real part of a CHEBFUN.
        f = real(f)
        
        % Restrict a CHEBFUN object to a subdomain.
        f = restrict(f, newDomain);

        % The roots of the CHEBFUN F.
        r = roots(f, varargin);
        
        % Round a CHEBFUN towards nearest integer.
        g = round(f)

        % Plot a CHEBFUN object on a log-linear scale:
        h = semilogx(f, varargin);

        % Plot a CHEBFUN object on a linear-log scale:
        h = semilogy(f, varargin);

        % Simplify the representation of a CHEBFUN obect.
        f = simplify(f, tol);

        % Size of a CHEBFUN object.
        [s1, s2] = size(f, dim);
        
        % Retrieve and modify preferences for this class.
        varargout = subsref(f, index);

        % Retrieve and modify preferences for this class.
        varargout = subsasgn(f, varargin);
        
        % CHEBFUN multiplication.
        f = times(f, g, varargin)
        
        % Transpose a CHEBFUN.
        f = transpose(f)
        
        % Adjust nearby common break points in domains of CHEBFUN objects.
        [f, g, newBreaksLocF, newBreaksLocG] = tweakDomain(f, g, tol)

        % Unary minus of a CHEBFUN.
        f = uminus(f)

        % Unary plus of a CHEBFUN.
        f = uplus(f)
        
        % Vertical scale of a CHEBFUN object.
        out = vscale(f);
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                (Private) Methods implemented in this m-file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function op = str2op(op)
% This is here as it's a clean function with no other variables hanging
% around in the scope.
depVar = symvar(op);
if ( numel(depVar) ~= 1 )
    error('CHEBFUN:STR2OP:indepvars', ...
        'Incorrect number of independent variables in string input.');
end
op = eval(['@(' depVar{:} ')', op]);
end

function [op, domain, pref] = parseInputs(op, domain, varargin)
% Parse inputs.

    args = varargin;
    if ( nargin == 1 )
        % chebfun(op)
        pref = chebpref();
        domain = pref.domain;
    elseif ( isstruct(domain) || isa(domain, 'chebpref') )
        % chebfun(op, pref)
        pref = chebpref(domain);
        domain = pref.domain;
    elseif ( ~isnumeric(domain) || (length(domain) == 1) )
        % chebfun(op, prop1, val1, ...)
        pref = chebpref();
        args = [domain, args];
        domain = pref.domain;
    elseif ( nargin < 3 )
        % chebfun(op, domain)
        pref = chebpref();
    elseif ( isstruct(varargin{1}) || isa(varargin{1}, 'chebpref') )
        % chebfun(op, domain, pref)
        pref = chebpref(args{1});
        args(1) = [];
    else
        % chebfun(op, domain, prop1, val1, ...)
        pref = chebpref();
    end

    % Take the default domain if an empty one was given:
    if ( isempty(domain) )
        domain = pref.domain;
    end

    vectorize = false;
    % Obtain additional preferences:
    while ( ~isempty(args) )
        if ( any(strcmpi(args{1}, {'chebpts1', 'chebpts2', 'equi'})) )
            % Determine tech for sampled values:
            if ( strcmpi(args{1}, 'chebpts1') )
                pref.tech = 'chebtech1';
            elseif ( strcmpi(args{1}, 'chebpts2') )
                pref.tech = 'chebtech2';
            elseif ( strcmpi(args{1}, 'equi') )
                pref.tech = 'funqui';
            end
            args(1) = [];
        elseif ( strcmpi(args{1}, 'vectorize') || ...
                 strcmpi(args{1}, 'vectorise') )
            % Vectorize flag for function_handles.
            vectorize = true;
            args(1) = [];
        elseif ( strcmpi(args{1}, 'coeffs') && isnumeric(op) )
            % Hack to support construction from coefficients.
            op = {{[], op}};
            args(1) = [];
        elseif ( isnumeric(args{1}) )
            % g = chebfun(@(x) f(x), N)
            pref.techPrefs.exactLength = args{1};
            args(1) = [];
        elseif ( strcmpi(args{1}, 'splitting') )
            % Translate "splitting" --> "enableBreakpointDetection".
            pref.enableBreakpointDetection = strcmpi(args{2}, 'on');
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'blowup') )
            % Translate "blowup" --> "enableSingularityDetection".
            pref.enableSingularityDetection = strcmpi(args{2}, 'on');
            args(1:2) = [];
        else
            % Update these preferences:
            pref.(args{1}) = args{2};
            args(1:2) = [];
        end
    end
    
    % [TODO]: Should we do the following for all elements in a cell input?
    % Convert string input to function_handle:
    if ( ischar(op) )
        op = str2op(op);
    end
    if ( isa(op, 'function_handle') && vectorize )
        % [TODO]: Should we reinstate VECTORCHECK()?
        op = vec(op);
    end
    
    if ( isa(op, 'function_handle') && strcmp(pref.tech, 'funqui') )
        if ( isfield(pref.techPrefs, 'exactLength') && ...
             ~isnan(pref.techPrefs.exactLength) )
            x = linspace(domain(1), domain(end), pref.techPrefs.exactLength).';
            op = feval(op, x);
            pref.techPrefs.exactLength = NaN;
        end
    end

end

function g = vec(f)
%VEC  Vectorize a function or string expression.
%   VEC(F), if F is a function handle or anonymous function, returns a function
%   that returns vector outputs for vector inputs by wrapping F inside a loop.
g = @loopwrapper;
    % Nested function:
    function v = loopwrapper(x)
        v = zeros(size(x));
        for j = 1:numel(v)
            v(j) = f(x(j));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% JUNK
% [TODO]: Reinstate or delete this.
%
% CHEBFUN(F,'extrapolate','on') prevents the constructor from evaluating
% the function F at the endpoints of the domain. This may also be achieved
% with CHEBFUN(F,'chebkind','1st','resampling','on') (which uses Chebyshev
% points of the 1st kind during the construction process), although this
% functionality is still experimental.
%
% CHEBFUN(F,...,'map',{MAPNAME,MAPPARS}) allows the use of mapped Chebyshev
% expansions. See help chebfun/maps for more information.
%
% CHEBFUN(CHEBS,ENDS,NP) specifies the number NP(i) of
% Chebyshev points for the construction of the function in CHEBS{i}.
%
% G = CHEBFUN(...) returns an object G of type chebfun.  A chebfun consists of a
% vector of 'funs', a vector 'domain' of length k+1 defining the intervals where
% the funs apply, and a matrix 'impulses' containing information about possible
% delta functions at the breakpoints between funs.

