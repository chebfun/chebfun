classdef trigtech < smoothfun % (Abstract)
%TRIGTECH   Approximate smooth periodic functions on [-1,1] with trigonometric 
%           interpolants.
%
%   Class for approximating smooth periodic functions on the interval [-1,1]
%   using function values at equally spaced points on [-1,1).
%
% Constructor inputs:
%   TRIGTECH(OP) constructs a TRIGTECH object from the function handle OP. OP
%   should be vectorized (i.e., accept a vector input) and ouput a vector of the
%   same length. TRIGTECH objects allow for array-valued construction (i.e., of
%   array-valued function), in which case OP should accept a vector of length N
%   and return a matrix of size NxM, where M is number of columns of the multi
%   -valued function.
%
%   TRIGTECH(OP, DATA) constructs a TRIGTECH using the additional data
%   supplied in the DATA structure.  Fields currently recognized are:
%     DATA.VSCALE    (Default:  0)
%     DATA.HSCALE    (Default:  1)
%         The constructor builds a TRIGTECH with 'happiness' (see
%         HAPPINESSCHECK.m) relative to the maximum of the given vertical scale
%         DATA.VSCALE and the (column-wise) infinity norm of the sampled
%         function values of OP, and the fixed horizontal scale DATA.HSCALE.
%   If any fields in DATA are empty or not supplied, or if DATA itself is empty
%   or not supplied, appropriate default values are set.
%
%   TRIGTECH(OP, DATA, PREF) overrides the default behavior with that given by
%   the preference structure PREF.
%
%   TRIGTECH(VALUES, ...) returns a TRIGTECH object which interpolates the
%   values in the columns of VALUES at equally spaced points and
%   TRIGTECH({VALUES, COEFFS}, ... ) uses the trig-coefficients passed in
%   COEFFS rather than computing them. If COEFFS are passed, the resulting
%   TRIGTECH is always deemed 'happy'.
%
% Examples: % Basic construction: f = trigtech(@(x) exp(sin(pi*x)))
%
%   % Construction with preferences:
%   p.sampleTest = 0; % See TRIGTECH.TECHPREF for details
%   f = trigtech(@(x) sin(x), [], [], p)
%
%   % Array-valued construction:
%   f = trigtech(@(x) tanh([sin(pi*x), cos(pi*x), cos(pi*sin(pi*x))]))
%
% See also TRIGTECH.TECHPREF, TRIGPTS, HAPPINESSCHECK, REFINE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIGTECH Class Description:
%
% The TRIGTECH class is for representations of smooth periodic functions on the
% interval [-1,1] via interpolated function values at equally spaced points
% using trigonometric Fourier series.
%
% The vertical scale VSCALE is used to enforce scale invariance in TRIGTECH
% construction and subsequent operations. For example, that
%   trigtech(@(x) 2^300*f(x)) = 2^300*trigtech(@(x) f(x)).
%
% VSCALE may be optionally passed to the constructor (if not, it defaults to 0),
% and during construction it is updated to be the maximum magnitude of the
% sampled function values. Similarly the horizontal scale HSCALE is used to
% enforce scale invariance when the input OP has been implicitly mapped from a
% domain other than [-1 1] before being passed to the TRIGTECH constructor.
%
% EPSLEVEL is the happiness level to which the TRIGTECH was constructed (See
% HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate of
% subsequent operations, both relative to VSCALE. Therefore EPSLEVEL could be
% regarded as the number of correct digits in the sampled value that created
% VSCALE.
%
% Here is a rough guide to how scale and accuracy information is propagated in
% subsequent operations after construction:
%   h = f + c:
%     h.vscale = max(abs(f.values), [], 1);
%     h.epslevel = (f.epslevel*f.vscale + eps(c)) / h.vscale;
%
%   h = f * c:
%     h.vscale = abs(c)*f.vscale;
%     h.epslevel = f.epslevel + eps(c)/c;
%
%   h = f + g:
%     h.vscale = max(abs(h.values), [], 1);
%     h.epslevel = (f.epslevel*f.vscale + g.epslevel*g.vscale) / h.vscale
%
%   h = f .* g:
%     h.vscale = max(abs(h.values), [], 1);
%     h.epslevel = (f.epslevel + g.epslevel) * (f.vscale*g.vscale)/h.vscale
%
%   h = diff(f):
%     h.vscale = max(abs(h.values), [], 1);
%     % [TODO]: Figure this out rigourously.
%     h.epslevel = n*log(n)*f.epslevel*f.vscale; % *(h.vscale/h.vscale)
%     % Note we don't divide by h.vscale here as we must also multiply by it.
%
%   h = cumsum(f):
%     h.vscale = max(abs(h.values), [], 1);
%     h.epslevel = happinessCheck(h);
%
% If the input operator OP in a call to TRIGTECH evaluates to NaN or Inf at
% any of the sample points used by the constructor, then an error is thrown.
%
% The TRIGTECH class supports the representation of array-valued functions (for
% example, f = trigtech(@(x) [sin(pi*x), cos(pi*x)])). In such cases, the values
% and coefficients are stored in a matrix (column-wise), and as such each
% component of the array-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All TRIGTECH methods should
% accept such array-valued forms. Note that this representation is distinct from
% an array of TRIGTECH objects, for which there is little to no support.
%
% Class diagram: [<<SMOOTHFUN>>] <-- [TRIGTECH]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )

        % Values of TRIGTECH at equally spaced points from [-1,1). For
        % array-valued TRIGTECH objects, each column represents the interpolated
        % values of a single function.
        values % (nxm double)

        % Coefficients are represented for the complex exponential form of the
        % interpolant. The coefficients are stored in descending order so that
        % c_{(N-1)/2} is the first entry and c_{-(N-1)/2} is the last. For
        % array-valued TRIGTECH objects, each column represents the coefficients
        % of a single function.
        coeffs % (nxm double)

        % Vertical scale of the TRIGTECH. This is a row vector storing the
        % magnitude of the largest entry in each column of VALUES. It is
        % convenient to store this as a property.
        vscale = 0 % (1xm double >= 0)

        % Horizontal scale of the TRIGTECH. Although TRIGTECH objects have in
        % principle no notion of horizontal scale invariance (since they always
        % live on [-1,1)), the input OP may have been implicitly mapped. HSCALE
        % is then used to enforce horizontal scale invariance in construction
        % and other subsequent operations that require it. It defaults to 1 and
        % is never updated.
        hscale = 1 % (scalar > 0)

        % Boolean value designating whether the TRIGTECH is 'happy' or not.
        % See HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)

        % Happiness level to which the TRIGTECH was constructed (See
        % HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate
        % of subsequent operations (See TRIGTECH class documentation for
        % details).
        epslevel % (double >= 0)
        
        % Boolean value designating whether the TRIGTECH represents a
        % real-valued function. This allows us to always return a real result
        % for things like evaluating a TRIGTECH.
        isReal % (logical)
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function obj = trigtech(op, data, pref)
            %Constructor for the TRIGTECH class.

            % Parse inputs.
            if ( (nargin == 0) || isempty(op) )
                % Return an empty TRIGTECH on null input:
                return
            end

            if ( (nargin < 2) || isempty(data) )
                    data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = trigtech.techPref();
            else
                pref = trigtech.techPref(pref);
            end

            data = parseDataInputs(data, pref);

            % Force nonadaptive construction if PREF.FIXEDLENGTH is numeric:
            if ( ~isempty(pref.fixedLength) && ~isnan(pref.fixedLength) )
                % Evaluate op on the equi-spaced grid of given size:
                vals = feval(op, trigtech.trigpts(pref.fixedLength));
                vals(1,:) = 0.5*(vals(1,:) + feval(op, 1));
                op = vals;
            end

            % Actual construction takes place here:
            obj = populate(obj, op, data.vscale, data.hscale, pref);
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % Absolute value of a TRIGTECH. (f should have no zeros in its domain)
        f = abs(f, pref)

        % TRIGTECH logical AND.
        h = and(f, g)

        % True if any element of a TRIGTECH is a nonzero number, ignoring NaN.
        a = any(f, dim)

        % Convert an array of TRIGTECH objects into an array-valued TRIGTECH.
        f = cell2mat(f)

        % Circular convolution of two trigtech objects.
        f = circconv(f, g)
        
        % Circular shift of a trigtech object by a: f -> g(x-a).
        g = circshift(f, a)

        % Plot (semilogy) the trigonometric coefficients of a TRIGTECH object.
        [h1, h2] = coeffsplot(f, varargin)

        % Check the happiness of a TRIGTECH. (Classic definition).
        [ishappy, epslevel, cutoff] = classicCheck(f, values, pref)

        % Compose two TRIGTECH objects or a TRIGTECH with a function handle:
        h = compose(f, op, g, data, pref)

        % Complex conjugate of a TRIGTECH.
        f = conj(f)
        
        % TRIGTECH objects are not transposable.
        f = ctranspose(f)

        % Indefinite integral of a TRIGTECH.
        f = cumsum(f, dim)

        % Derivative of a TRIGTECH.
        f = diff(f, k, dim)
        
        % Extract information for DISPLAY.
        info = dispData(f)
        
        % Extract columns of an array-valued TRIGTECH object.
        f = extractColumns(f, columnIndex)

        % Evaluate a TRIGTECH.
        y = feval(f, x)
        
        % Flip columns of an array-valued TRIGTECH object.
        f = fliplr(f)
        
        % Flip/reverse a TRIGTECH object.
        f = flipud(f)

        % Plot (semilogy) the trigonometric coefficients of a TRIGTECH object.
        varargout = plotcoeffs(f, varargin)

        % Get method:
        val = get(f, prop);

        % Happiness test for a TRIGTECH
        [ishappy, epslevel, cutoff] = happinessCheck(f, op, values, pref)

        % Imaginary part of a TRIGTECH.
        f = imag(f)

        % Compute the inner product of two TRIGTECH objects.
        out = innerProduct(f, g)

        % True for an empty TRIGTECH.
        out = isempty(f)

        % Test if TRIGTECH objects are equal.
        out = isequal(f, g)

        % Test if a TRIGTECH is bounded.
        out = isfinite(f)

        % Test if a TRIGTECH is unbounded.
        out = isinf(f)

        % Test if a TRIGTECH has any NaN values.
        out = isnan(f)

        function out = isPeriodicTech(f)
        %ISPERIODICTECH   True for TRIGTECH.
            out = 1;
        end
        
        % True for real TRIGTECH.
        out = isreal(f)
        
        % True for zero TRIGTECH objects
        out = iszero(f)
        
        % Cannot convert TRIGTECH coefficients to legendre coefficients
        c_leg = legcoeffs(f, n)
        
        % Length of a TRIGTECH.
        len = length(f)

        % Convert an array-valued TRIGTECH into an ARRAY of TRIGTECH objects.
        g = mat2cell(f, M, N)

        % Global maximum of a TRIGTECH on [-1,1].
        [maxVal, maxPos] = max(f)

        % Global minimum of a TRIGTECH on [-1,1].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [-1,1].
        [vals, pos] = minandmax(f)

        % Subtraction of two TRIGTECH objects.
        f = minus(f, g)

        % Left matrix divide for TRIGTECH objects.
        X = mldivide(A, B)

        % Right matrix divide for a TRIGTECH.
        X = mrdivide(B, A)

        % Multiplication of TRIGTECH objects.
        f = mtimes(f, c)

        % TRIGTECH logical OR.
        h = or(f, g)

        % Basic linear plot for TRIGTECH objects.
        varargout = plot(f, varargin)
        
        % 3-D plot for TRIGTECH objects.
        varargout = plot3(f, g, h, varargin)
        
        % Obtain data used for plotting a TRIGTECH object:
        data = plotData(f, g, h)

        % Addition of two TRIGTECH objects.
        f = plus(f, g)

        % Return the points used by a TRIGTECH.
        out = points(f)

        % Complex polynomial coefficients of a TRIGTECH.
        out = poly(f)

        % Populate a TRIGTECH class with values.
        f = populate(f, op, vscale, hscale, pref)
        
        % Power function of a TRIGTECH.
        f = power(f, b)

        % Adjust the number of points used in a TRIGTECH.
        f = prolong(f, n)

        % QR factorisation of an array-valued TRIGTECH.
        [f, R, E] = qr(f, flag, methodFlag)

        % Right array divide for a TRIGTECH.
        f = rdivide(f, c, pref)

        % Real part of a TRIGTECH.
        f = real(f)

        % Restrict a TRIGTECH to a subinterval.
        f = restrict(f, s)

        % Roots of a TRIGTECH in the interval [-1,1].
        out = roots(f, varargin)
        
        % Test an evaluation of the input OP against a TRIGTECH approx.
        pass = sampleTest(op, values, f)
        
        % Signum of a TRIGTECH. (f should have no zeros in its domain)
        f = sign(f, pref)

        % Trim trailing Chebyshev coefficients of a TRIGTECH object.
        f = simplify(f, pref, force)

        % Size of a TRIGTECH.
        [siz1, siz2] = size(f, varargin)

        % Definite integral of a TRIGTECH on the interval [-1,1].
        out = sum(f, dim)

        % TRIGTECH multiplication.
        f = times(f, g, varargin)
        
        % TRIGTECH obects are not transposable.
        f = transpose(f)
        
        % Coefficients for degree N least square approximation of F:
        c = truncCoeffs(f, n)

        % Unary minus of a TRIGTECH.
        f = uminus(f)

        % Unary plus of a TRIGTECH.
        f = uplus(f)   
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static = true )
        
        % Aliasing:
        coeffs = alias(coeffs, m)
        
        % Differentiation matrix in Fourier basis.
        D = diffmat(n, p)
        
        % Compute trigonometric points (x) and optionally quadrature (w)
        % and barycentric (v) weights:
        [x, w] = trigpts(n);
        
        % Convert coefficients to values:
        values = coeffs2vals(coeffs);
        
        % Make a TRIGTECH (constructor shortcut):
        f = make(varargin);
        
        % Compute trigonometric quadrature weights (trapezoidal rule):
        w = quadwts(n)
        
        % Refinement function for TRIGTECH construction (evaluates OP on grid):
        [values, points, giveUp] = refine(op, values, pref)
        
        % Retrieve and modify preferences for this class.
        p = techPref(q)

        % Convert values to coefficients:
        coeffs = vals2coeffs(values)

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( ~isfield(data, 'vscale') || isempty(data.vscale) )
    data.vscale = 0;
end

if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    data.hscale = 1;
end

end
