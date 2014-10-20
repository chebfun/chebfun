classdef chebtech < smoothfun % (Abstract)
%CHEBTECH   Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Abstract class for approximating smooth functions on the interval [-1,1]
%   using function values at Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% See also CHEBTECH1, CHEBTECH2, CHEBTECH.TECHPREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBTECH Class Description:
%
% The CHEBTECH class is an abstract class for representations of smooth
% functions on the interval [-1,1] via interpolated function values at Chebyshev
% points and coefficients of the corresponding first-kind Chebyshev series
% expansion.
%
% There are two concrete realizations of the CHEBTECH class--CHEBTECH1 and
% CHEBTECH2--which interpolate on Chebyshev grids of the 1st and 2nd kind,
% respectively.  Note that although they use different Chebyshev grids in
% 'value' space, their coefficients are always from an expansion in first-kind
% Chebyshev polynomials (i.e., those usually denoted by $T_k(x)$).
%
% The vertical scale VSCALE is used to enforce scale invariance in CHEBTECH
% construction and subsequent operations. For example, that
%   chebtech1(@(x) 2^300*f(x)) = 2^300*chebtech1(@(x) f(x)).
%
% VSCALE may be optionally passed to the constructor (if not, it defaults to 0),
% and during construction it is updated to be the maximum magnitude of the
% sampled function values. Similarly the horizontal scale HSCALE is used to
% enforce scale invariance when the input OP has been implicitly mapped from a
% domain other than [-1 1] before being passed to a CHEBTECH constructor.
%
% EPSLEVEL is the happiness level to which the CHEBTECH was constructed (See
% HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate of
% subsequent operations, both relative to VSCALE. Therefore EPSLEVEL could be
% regarded as the number of correct digits in the sampled value that created
% VSCALE.
%
% Here is a rough guide to how scale and accuracy information is propagated in
% subsequent operations after construction:
%   h = f + c:
%     h.vscale = getvscl(h);
%     h.epslevel = (f.epslevel*f.vscale + eps(c)) / h.vscale;
%
%   h = f * c:
%     h.vscale = getvscl(h) = abs(c)*f.vscale;
%     h.epslevel = f.epslevel + eps(c)/c;
%
%   h = f + g:
%     h.vscale = getvscl(h);
%     h.epslevel = (f.epslevel*f.vscale + g.epslevel*g.vscale) / h.vscale
%
%   h = f .* g:
%     h.vscale = getvscl(h);
%     h.epslevel = (f.epslevel + g.epslevel) * (f.vscale*g.vscale)/h.vscale
%
%   h = diff(f):
%     h.vscale = getvscl(h);
%     % [TODO]: Figure this out rigourously.
%     h.epslevel = n*log(n)*f.epslevel*f.vscale; % *(h.vscale/h.vscale)
%     % Note we don't divide by h.vscale here as we must also multiply by it.
%
%   h = cumsum(f):
%     h.vscale = getvscl(h);
%     % [TODO]: Figure this out rigourously.
%     h.epslevel = 2*f.epslevel*f.vscale/h.vscale
%
% If the input operator OP in a call to a concrete CHEBTECH constructor, say,
% CHEBTECH1(OP), evaluates to NaN or Inf at any of the sample points used by the
% constructor, then a suitable replacement is found by extrapolating (globally)
% from the numeric values (see EXTRAPOLATE.M). If the EXTRAPOLATE preference is
% set to TRUE (See CHEBTECH.TECHPREF), then the endpoint values -1 and +1 are
% always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The CHEBTECH classes support the representation of array-valued functions (for
% example, f = chebtech1(@(x) [sin(x), cos(x)])). In such cases, the values and
% coefficients are stored in a matrix (column-wise), and as such each component
% of the array-valued function is truncated to the same length, even if the
% demands of 'happiness' imply that one of the components could be truncated to
% a shorter length than the others. All CHEBTECH methods should accept such
% array-valued forms. Note that this representation is distinct from an array of
% CHEBTECH objects, for which there is little to no support.
%
% Class diagram: [<<SMOOTHFUN>>] <-- [<<CHEBTECH>>] <-- [CHEBTECH1]
%                                                   <-- [CHEBTECH2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )

        % Coefficients in 1st-kind Chebyshev series expansion of the CHEBTECH on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last. For array-valued CHEBTECH
        % objects, each column represents the coefficients of a single function.
        coeffs % (nxm double)

        % Vertical scale of the CHEBTECH. This is a row vector storing the
        % magnitude of the largest entry in each column of VALUES. It is
        % convenient to store this as a property.
        vscale = 0 % (1xm double >= 0)

        % Horizontal scale of the CHEBTECH. Although CHEBTECH objects have in
        % principle no notion of horizontal scale invariance (since they always
        % live on [-1,1]), the input OP may have been implicitly mapped. HSCALE
        % is then used to enforce horizontal scale invariance in construction
        % and other subsequent operations that require it. It defaults to 1 and
        % is never updated.
        hscale = 1 % (scalar > 0)

        % Boolean value designating whether the CHEBTECH is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)

        % Happiness level to which the CHEBTECH was constructed (See
        % HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate
        % of subsequent operations (See CHEBTECH class documentation for
        % details).
        epslevel % (double >= 0)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Abstract = true )

        % Compose method. (Not implemented here as refinement is defined also).
        h = compose(f, op, g, data, pref)

        % Get method.
        val = get(f, prop);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true, Abstract = true )
        
        % Alias Chebyshev coefficients.
        coeffs = alias(coeffs, m)
        
        % Angles of Chebyshev points. (i.e., acos(chebpts(n))
        t = angles(n)

        % Compute Chebyshev barycentric weights.
        w = barywts(n)

        % Compute Chebyshev points (x) and optionally quadrature (w) and
        % barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Convert coefficients to values.
        values = coeffs2vals(coeffs)

        % Make a CHEBTECH. (Constructor shortcut)
        f = make(varargin);

        % Refinement function for CHEBTECH construction. (Evaluates OP on grid)
        [values, points, giveUp] = refine(op, values, pref)

        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Convert values to coefficients.
        coeffs = vals2coeffs(values)

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % Absolute value of a CHEBTECH. (f should have no zeros in its domain)
        f = abs(f, pref)

        % CHEBTECH logical AND.
        h = and(f, g)

        % True if any element of a CHEBTECH is a nonzero number, ignoring NaN.
        a = any(f, dim)

        % Convert an array of CHEBTECH objects into an array-valued CHEBTECH.
        f = cell2mat(f)

        % Plot (semilogy) the Chebyshev coefficients of a CHEBTECH object.
        [h1, h2] = plotcoeffs(f, varargin)

        % Check the happiness of a CHEBTECH. (Classic definition).
        [ishappy, epslevel, cutoff] = classicCheck(f, values, pref)

        % Complex conjugate of a CHEBTECH.
        f = conj(f)
        
        % CHEBTECH objects are not transposable.
        f = ctranspose(f)

        % Indefinite integral of a CHEBTECH.
        f = cumsum(f, dim)

        % Derivative of a CHEBTECH.
        f = diff(f, k, dim)

        % Extract information for DISPLAY.
        info = dispData(f)
        
        % Extract columns of an array-valued CHEBTECH object.
        f = extractColumns(f, columnIndex)

        % Extract roots at the boundary points -1 and 1.
        [f, rootsLeft, rootsRight] = extractBoundaryRoots(f, numRoots)

        % Extrapolate (for NaNs / Infs).
        [values, maskNaN, maskInf] = extrapolate(f, values)

        % Evaluate a CHEBTECH.
        y = feval(f, x)
        
        % Round a CHEBTECH towards zero.
        g = fix(f);
        
        % Round a CHEBTECH towards minus infinity.
        g = floor(f);

        % Flip columns of an array-valued CHEBTECH object.
        f = fliplr(f)
        
        % Flip/reverse a CHEBTECH object.
        f = flipud(f)

        % Happiness test for a CHEBTECH
        [ishappy, epslevel, cutoff] = happinessCheck(f, op, values, pref)

        % Imaginary part of a CHEBTECH.
        f = imag(f)

        % Compute the inner product of two CHEBTECH objects.
        out = innerProduct(f, g)

        % Test if a CHEBTECH decays faster than a single root at endpoints.
        out = isdecay(f)

        % True for an empty CHEBTECH.
        out = isempty(f)

        % Test if CHEBTECH objects are equal.
        out = isequal(f, g)

        % Test if a CHEBTECH is bounded.
        out = isfinite(f)

        % Test if a CHEBTECH is unbounded.
        out = isinf(f)

        % Test if a CHEBTECH has any NaN values.
        out = isnan(f)

        % True for real CHEBTECH.
        out = isreal(f)
        
        % True for zero CHEBTECH objects
        out = iszero(f)
        
        % Return Legendre coefficients of a CHEBTECH.
        c_leg = legcoeffs(f, n)
        
        % Length of a CHEBTECH.
        len = length(f)

        % CHEBTECH logical.
        f = logical(f)

        % A 'loose' (i.e., not too strict) check for happiness.
        [ishappy, epslevel, cutoff] = looseCheck(f, values, pref)

        % Convert an array-valued CHEBTECH into an ARRAY of CHEBTECH objects.
        g = mat2cell(f, M, N)

        % Global maximum of a CHEBTECH on [-1,1].
        [maxVal, maxPos] = max(f)

        % Global minimum of a CHEBTECH on [-1,1].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [-1,1].
        [vals, pos] = minandmax(f)

        % Subtraction of two CHEBTECH objects.
        f = minus(f, g)

        % Left matrix divide for CHEBTECH objects.
        X = mldivide(A, B)

        % Right matrix divide for a CHEBTECH.
        X = mrdivide(B, A)

        % Multiplication of CHEBTECH objects.
        f = mtimes(f, c)

        % CHEBTECH logical NOT.
        f = not(f)

        % CHEBTECH logical OR.
        h = or(f, g)

        % Basic linear plot for CHEBTECH objects.
        varargout = plot(f, varargin)
        
        % 3-D plot for CHEBTECH objects.
        varargout = plot3(f, g, h, varargin)
        
        % Obtain data used for plotting a CHEBTECH object:
        data = plotData(f, g, h)

        % Addition of two CHEBTECH objects.
        f = plus(f, g)

        % Return the points used by a CHEBTECH.
        out = points(f)

        % Polynomial coefficients of a CHEBTECH.
        out = poly(f)

        % Populate a CHEBTECH class with values.
        [f, values] = populate(f, op, vscale, hscale, pref)
        
        % Power function of a CHEBTECH.
        f = power(f, b)

        % Adjust the number of points used in a CHEBTECH.
        f = prolong(f, n)

        % QR factorisation of an array-valued CHEBTECH.
        [f, R, E] = qr(f, flag, methodFlag)

        % Right array divide for a CHEBTECH.
        f = rdivide(f, c, pref)

        % Real part of a CHEBTECH. REQUIRED BY THIS CLASS:
        f = real(f)

        % Restrict a CHEBTECH to a subinterval.
        f = restrict(f, s)

        % Roots of a CHEBTECH in the interval [-1,1].
        out = roots(f, varargin)
        
        % Round a CHEBTECH towards nearest integer.
        g = round(f)

        % Test an evaluation of the input OP against a CHEBTECH approx.
        pass = sampleTest(op, values, f)
        
        % Signum of a CHEBTECH. (f should have no zeros in its domain)
        f = sign(f, pref)

        % Trim trailing Chebyshev coefficients of a CHEBTECH object.
        f = simplify(f, pref, force)

        % Size of a CHEBTECH.
        [siz1, siz2] = size(f, varargin)

        % Strict happiness check.
        [ishappy, epslevel, cutoff] = strictCheck(f, values, pref)

        % Definite integral of a CHEBTECH on the interval [-1,1].
        out = sum(f, dim)

        % CHEBTECH multiplication.
        f = times(f, g, varargin)
         
        % CHEBTECH obects are not transposable.
        f = transpose(f)

        % Unary minus of a CHEBTECH.
        f = uminus(f)

        % Unary plus of a CHEBTECH.
        f = uplus(f)

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS (IMPLEMENTED BY THIS CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )

        % Clenshaw's algorithm for evaluating a Chebyshev polynomial.
        out = clenshaw(x, coeffs)

        % Retrieve and modify preferences for this class.
        p = techPref(q)

    end

end
