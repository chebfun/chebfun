classdef funcheb2 < funcheb
%FUNCHEB2   Approximate smooth functions on [-1,1] with Chebyshev interpolants. 
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values at 2nd-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   FUNCHEB2(OP) constructs a FUNCHEB2 object from the function handle OP. OP
%   should be vectorised (i.e., accept a vector input) and ouput a vector of the
%   same length. FUNCHEB2 objects allow for vectorised construction (i.e., of a
%   multi-valued function), in which case OP should accept a column vector of
%   length N and return a matrix of size NxM.
%
%   FUNCHEB2(OP, VSCALE) constructs a FUNCHEB2 with 'happiness' (see
%   HAPPINESSCHECK.m) relative to the maximum of the given vertical scale VSCALE
%   and the infinity norm of the sampled function values of OP. If not given,
%   the VSCALE defaults to 0 initially.
%
%   FUNCHEB2(OP, VSCALE, PREF) overrides the default behavior with that given by
%   the preference structure PREF. The constructor will also accept inputs of
%   the form FUNCHEB2(OP, PREF), but this usage is not advised. Similarly, one
%   can pass FUNCHEB2(OP, VSCALE, EPS), which is equivalent to the call
%   FUNCHEB2(OP, VSCALE, FUNCHEB2.PREF('eps',EPS)).
%
%   FUNCHEB2(VALUES) or FUNCHEB2(VALUES, VSCALE, PREF) returns a FUNCHEB2 object
%   which interpolates the values in the columns of VALUES at 2nd-kind Chebyshev
%   points.
%
% Examples:
%   % Basic construction:
%   f = funcheb2(@(x) sin(x))
%
%   % Construction with preferences:
%   p = funcheb2.pref('sampletest', 0); % See help('funcheb2.pref') for details
%   f = funcheb2(@(x) sin(x), p)
%
%   % Vector-valued construction:
%   f = funcheb2(@(x) [sin(x), cos(x), exp(x)])
%
% See also FUNCHEB2.pref, FUNCHEB2.chebpts, FUNCHEB2.happinesscheck, FUNCHEB2.refine.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB2 Class Description:
%
% The FUNCHEB2 class represents smooth functions on the interval [-1,1] using
% function values at 2nd-kind Chebyshev points and coefficients of the
% corresponding first-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 2nd-kind grid (see REFINE.m) until the
% representation is deemed 'happy' (see HAPPINESSCHECK.m). The resulting object
% can be used to evaluate and operate on the input function.
%
% The vertical scale VSCALE is used to enforce scale invariance in FUNCHEB2
% construction and subsequent operations. For example, that 
%          funcheb2(@(x) 2^300*f(x)) = 2^300*funcheb2(@(x) f(x)). 
% VSCALE may be optionally passed during to the constructor (if not, it defaults
% to 0), and during construction it is updated to be the maximum magnitude of
% the sampled function values.
%
% FUNCHEB2 objects have no notion of horizontal scale invariance (since they
% always live on [-1,1]). However, if the input OP has been implicitly mapped
% one can enforce construction relative to some horizontal scale HSCALE by using
% the preference FUNCHEB2.pref('hscale', HSCALE).
%
% If the input operator OP evaluates to NaN or Inf at any of the sample points
% used by the constructor, then a suitable replacement is found by extrapolating
% (globally) from the numeric values (see EXTRAPOLATE.M). If the preference
% funcheb2.pref('extrapolate', TRUE) is set, then the endpoint values -1 and +1
% are always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The FUNCHEB2 class supports the representation of vector-valued functions
% (for example, f = funcheb2(@(x) [sin(x), cos(x)])). In such cases, the values
% and coefficients are stored in a matrix (column-wise), and as such each
% component of the multi-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All FUNCHEB2 methods should
% accept such vectorised forms. Note that this representation is distinct from
% an array of funcheb2 objects, such as [funcheb2(@(x) sin(x), funcheb2(@(x)
% cos(x)], for which there is little to no support.
%
% Class diagram: [<<funcheb>>] <-- [FUNCHEB2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Properties of FUNCHEB2 objects.
    properties ( Access = public )
        
        % Values of FUNCHEB2 at 2nd-kind Chebyshev points. Values are stored in
        % order from left to right (i.e., values(1,:) = feval(FUNCHEB2, -1)).
        % For vectorised FUNCHEB2 objects, each column represents the
        % interpolated values of a single function.
        values % (nxm double)
        
        % Coefficients in 1st-kind Chebyshev series expansion of the FUNCHEB2 on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last. For vectorised FUNCHEB2 objects,
        % each column represents the coefficients of a single function.
        coeffs % (nxm double)
        
        % Vertical scale of the FUNCHEB2. Typically the is the magnitude of the
        % largest value sampled from the given operator OP during the
        % construction process. It is updated via subsequent FUNCHEB2
        % operations in a natural way.
        vscale % (1xm double >= 0)
        
        % Boolean value designating whether the FUNCHEB2 is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)
        
        % Happiness level to which the FUNCHEB2 was constructed, or a rough
        % accuracy estimate of subsequent operations. See HAPPINESSCHECK.m for
        % full documentation.
        epslevel % (double >= 0)
    end
    
    %% Methods implemented by this m-file.
    methods
        
        function obj = funcheb2(op, vscale, pref)
            %Constructor for the FUNCHEB2 class.
            
            % Return an empty FUNCHEB2 on null input:
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
                pref = funcheb2.pref(vscale);
                vscale = 0;
            elseif ( nargin < 3 )
                % Create:
                pref = funcheb2.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed.
                pref = funcheb2.pref('eps', pref);
            else
                % Merge:
                pref = funcheb2.pref(pref);
            end
            
            % Force nonadaptive construction if pref.funcheb2.n is numeric:
            if ( ~isempty(pref.funcheb2.n) && ~isnan(pref.funcheb2.n) )
                % Evaluate the op on the Chebyshev grid of given size:
                op = feval(op, funcheb2.chebpts(pref.funcheb2.n));
            end
            
            % Actual construction takes place here:
            if ( ~isnumeric(op) )
                % Adaptive contruction:
                [values, coeffs, vscale, ishappy, epslevel] = ...
                    funcheb2.constructor(op, vscale, pref); %#ok<*PROP>
            else
                % Nonadaptive contruction: (op is a vector of values)
                values = op;
                coeffs = funcheb2.chebpoly(values);
                vscale = max(abs(values));
                [ishappy, epslevel] = ...
                    funcheb2.happinessCheck(op, values, coeffs, vscale, pref);
            end
            
            % Assign to FUNCHEB2 object:
            obj.values = values;
            obj.coeffs = coeffs;
            obj.vscale = vscale;
            obj.ishappy = ishappy;
            obj.epslevel = epslevel;
            
        end
        
    end
    
    %% Static methods implemented by FUNCHEB2 class.
    % (This list is alphabetical)
    methods ( Static = true )
        
        % Alias Chebyshev coefficients.
        coeffs = alias(coeffs, m);
        
        % Evaluate a Chebyshev interpolant using barycentric formula.
        out = bary(x, values, kind)
        
        % Compute Chebyshev barycentric weights.
        w = barywts(n)
        
        % Convert values to coefficients.
        coeffs = chebpoly(values)
        
        % Convert coefficients to values.
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points of first kind (x) and optionally quadrature
        % (w) and barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % The constructor for the FUNCHEB2 class.
        [values, coeffs, vscale, ishappy, epslevel] = ...
            constructor(op, vscale, pref)
        
        % Evaluate a Chebyshev polynomial using Clenshaw's algorithm.
        out = clenshaw(x, coeffs)
        
        % Default happiness in FUNCHEB2 (from Chebfun v4).
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
        
        % Refinement function for FUNCHEB2 construction. (Evaluates OP on grid)
        [values, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Test a sample evaluation of an interpolant against op evaluation.
        pass = sampleTest(op, values, vscale, epslevel, pref)
        
        % Strict happiness check.
        [ishappy, epslevel, cotoff] = ...
            strictCheck(op, values, coeffs, vscale, pref)
        
        % Test the FUNCHEB2 class.
        pass = test(varargin);
        
    end
    
    %% Private static methods implemented by FUNCHEB2 class.
    methods ( Access = public, Static = true )
        
    end
    
end
