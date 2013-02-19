classdef funcheb1 < funcheb
%FUNCHEB1 Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Class for approximating smooth functions on the interval [-1,1]
%   using function values at 1st-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   FUNCHEB1(OP) constructs a FUNCHEB1 object from the function handle OP. OP
%   should be vectorized (i.e., accept a vector input) and ouput a vector of
%   the same length. FUNCHEB1 objects allow for vectorised construction (i.e.,
%   of multi-valued function), in which case OP should accept a vector of length
%   N and return a matrix of size NxM.
%
%   FUNCHEB1(OP, VSCALE) constructs a FUNCHEB1 with 'happiness' (see
%   funcheb1.HAPPINESSCHECK.m) relative to the maximum of the given vertical scale (VSCALE)
%   and the infinity norm of the sampled function values of OP. If not given,
%   the VSCALE defaults to 0 initially.
%
%   FUNCHEB1(OP, VSCALE, PREF) overrides the default behavior with that given by
%   the preference structure PREF. The constructor will also accept inputs of
%   the form FUNCHEB2(OP, PREF), but this usage is not advised. Similarly, one
%   can pass FUNCHEB2(OP, VSCALE, EPS), which is equivalent to the call
%   FUNCHEB1(OP, VSCALE, FUNCHEB1.PREF('eps',EPS)).
%
%   FUNCHEB1(VALUES) or FUNCHEB1(VALUES, VSCALE, PREF) returns a FUNCHEB2 object
%   which interpolates the values in the columns of VALUES at 1st-kind Chebyshev
%   points.
%
% Examples:
%   % Basic construction
%   f = funcheb1(@(x) sin(x))
%
%   % Construction with preferences
%   p = funcheb1.pref('sampletest', 0); % See help('funcheb1.pref') for details
%   f = funcheb1(@(x) sin(x), p)
%
%   % Vector-valued construction:
%   f = funcheb1(@(x) [sin(x), cos(x), exp(x)])
%
% See also FUNCHEB1.pref, FUNCHEB1.chebpts, FUNCHEB1.happinesscheck, FUNCHEB1.refine.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB1 Class Description:
%
% The FUNCHEB1 class represents smooth functions on the interval [-1,1] using
% function values at 1st-kind Chebyshev points and coefficients of the
% corresponding first-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 1st-kind grid (see funcheb1.REFINE.m) until the
% representation is deemed 'happy' (see funcheb1.HAPPINESSCHECK.m). The resulting
% object can be used to evaluate and operate on the input function.
%
% The vertical scale VSCALE is used to enforce scale invariance in FUNCHEB1
% construction and subsequent operations. For example, that 
%          funcheb1(@(x) 2^300*f(x)) = 2^300*funcheb1(@(x) f(x)). 
% VSCALE may be optionally passed during to the constructor (if not, it
% defaults to 0), and during construction it is updated to be the maximum
% magnitude of the sampled function values.
%
% FUNCHEB1 objects have no notion of horizontal scale invariance (since they
% always live on [-1,1]). However, if the input OP has been implicitly mapped
% one can enforce construction relative to some horizontal scale HSCALE by using
% the preference FUNCHEB1.pref('hscale', HSCALE).
%
% If the input operator OP evaluates to NaN or Inf at any of the sample points
% used by the constructor, then a suitable replacement is found by extrapolating
% (globally) from the numeric values (see EXTRAPOLATE.M). If the preference
% funcheb1.pref('extrapolate', TRUE) is set, then the endpoint values -1 and +1
% are always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The FUNCHEB1 class supports the representation of vector-valued functions
% (for example, f = funcheb1(@(x) [sin(x), cos(x)])). In such cases, the values
% and coefficients are stored in a matrix (column-wise), and as such each
% component of the multi-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All FUNCHEB1 methods should
% accept such vectorised forms. Note that this representation is distinct from
% an array of funcheb1 objects, such as [funcheb1(@(x) sin(x), funcheb1(@(x)
% cos(x)], for which there is little to no support.
%
% Class diagram: [funcheb] <-- [funcheb1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Properties of FUNCHEB1 objects.
    properties ( Access = public )
        
        % Values of FUNCHEB1 at 1st-kind Chebyshev points. Values are stored in
        % order from left to right (i.e., values(1,:) = feval(FUNCHEB1,-1)).
        % For vectorised FUNCHEB1 objects, each column represents the
        % interpolated values of a single function.
        values % (n x m double)
        
        % Coefficients in 1st-kind Chebyshev series expansion of the FUNCHEB1 on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last. For vectorised FUNCHEB1 objects,
        % each column represents the coefficients of a single function.
        coeffs % (n x m double)
        
        % Vertical scale of the FUNCHEB1. Typically the is the magnitude of the
        % largest value sampled from the given operator OP during the
        % construction process. It is updated via subsequent FUNCHEB1
        % operations in a natural way.
        vscale % (double >= 0)
        
        % Boolean value designating whether the FUNCHEB1 is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)
        
        % Happiness level to which the FUNCHEB1 was constructed, or a rough
        % accuracy estimate of subsequent operations. See HAPPINESSCHECK.m for
        % full documentation.
        epslevel % (double >= 0)
    end
    
    %% Methods implemented by this m-file.
    methods
        
        function obj = funcheb1(op, vscale, pref)
            % Constructor for the FUNCHEB1 class.
            
            % Return an empty FUNCHEB1 on null input:
            if ( nargin == 0 || isempty(op) )
                return
            end
           
            % Define vscale if none given:
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            
            % Obtain preferences:
            if ( nargin == 2 && isstruct(vscale) )
                % vscale was actually a preference.
                pref = funcheb1.pref(vscale);
                vscale = 0;
            elseif ( nargin < 3 )
                % Create
                pref = funcheb1.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed
                pref = funcheb1.pref('eps', pref);
            else
                % Merge
                pref = funcheb1.pref(pref);
            end
            
            % Force nonadaptive construction if pref.funcheb1.n is numeric:
            if ( ~isempty(pref.funcheb1.n) && ~isnan(pref.funcheb1.n) )
                % Evaluate the op on the Chebyshev grid of given size:
                op = feval(op, funcheb1.chebpts(pref.funcheb1.n));
            end
            
            % Actual construction takes place here:
            if ( ~isnumeric(op) )
                % Adaptive contruction:
                [values, coeffs, vscale, ishappy, epslevel] = ...
                    funcheb1.constructor(op, vscale, pref); %#ok<*PROP>
            else
                % Nonadaptive contruction:
                values = op;
                coeffs = funcheb1.chebpoly(values);
                vscale = norm(values(:), inf);
                [ishappy, epslevel] = ...
                    funcheb1.happinessCheck(op, values, coeffs, vscale, pref);
            end
            
            % Assign to FUNCHEB1 object.
            obj.values = values;
            obj.coeffs = coeffs;
            obj.vscale = vscale;
            obj.ishappy = ishappy;
            obj.epslevel = epslevel;
            
        end
        
    end
    
    %% Static methods implemented by FUNCHEB1 class.
    % (This list is alphabetical)
    methods ( Static = true )
        
        % Alias Chebyshev coefficients.
        coeffs = alias(coeffs, m);
        
        % Evaluate a Chebyshev interpolant using barycentric formula.
        out = bary(x, values, kind)
        
        % Compute Chebyshev barycentric weights.
        w = barywts(n)
        
        % Barycentric Interpolation matrix using Chebyshev points of 1st 
        % kind.
        B = barymat(y, x, w)

        % Convert values to coefficients.
        coeffs = chebpoly(values)
        
        % Convert coefficients to values.
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points of first kind (x) and optionally
        % quadrature (w) and barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % The constructor for the FUNCHEB1 class.
        [values, coeffs, vscale, ishappy, epslevel] = ...
            constructor(op, vscale, pref)
        
        % Evaluate a Chebyshev polynomial using Clenshaw's algorithm.
        out = clenshaw(x, coeffs)
        
        % Default happiness in FUNCHEB1 (from Chebfun v4).
        [cutoff, epslevel] = classicCheck(values, coeffs, vscale, pref)       
        
        % Extrapolate (for NaNs / endpoints).
        [values, maskNaN, maskInf] = extrapolate(values)
        
        % Happiness checker.
        [ishappy, epslevel, cutoff] = ...
            happinessCheck(op, values, coeffs, vscale, pref)
        
        % Loose happiness check.
        [ishappy, epslevel, cotoff] = ...
            looseCheck(op, values, coeffs, vscale, pref)        
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Refinement function for FUNCHEB1 construction. (Evaluates OP on grid)
        [values, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Test a sample evaluation of an interpolant against op evaluation.
        pass = sampleTest(op, values, vscale, epslevel, pref)
        
        % Strict happiness check.
        [ishappy, epslevel, cotoff] = ...
            strictCheck(op, values, coeffs, vscale, pref)
        
        % Test the FUNCHEB1 class.
        pass = test(varargin);
        
    end
    
    %% Private static methods implemented by FUNCHEB1 class.
    methods ( Access = private, Static = true )
        
    end
    
end
