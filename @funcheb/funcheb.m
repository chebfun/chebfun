classdef funcheb < smoothfun % (Abstract)
%FUNCHEB   Approximate smooth functions on [-1,1] with Chebyshev interpolants. 
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values Chebyshev points and coefficients of the corresponding
%   1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   FUNCHEB.CONSTRUCTOR(OP) constructs a FUNCHEB object from the function handle
%   OP. OP should be vectorised (i.e., accept a vector input) and ouput a vector
%   of the same length. FUNCHEB objects allow for vectorised construction (i.e.,
%   of multi-valued functions), in which case OP should accept a column vector
%   of length N and return a matrix of size NxM.
%
%   FUNCHEB.CONSTRUCTOR(OP, VSCALE, HSCALE) constructs a FUNCHEB with
%   'happiness' relative to the maximum of the given vertical scale VSCALE
%   (which is updated by the infinity norm of the sampled function values of OP
%   during construction), and the fixed horizontal scale HSCALE. If not given,
%   the VSCALE defaults to 0 initially, and HSCALE defaults to 1.
%
%   FUNCHEB.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) overrides the default behavior
%   with that given by the preference structure PREF. The constructor will also
%   accept inputs of the form FUNCHEB(OP, PREF), but this usage is not advised.
%   Similarly, one can pass FUNCHEB(OP, VSCALE, HSCALE< EPS), which is
%   equivalent to the call FUNCHEB(OP, VSCALE, HSCALE, FUNCHEB.PREF('eps',EPS)).
%
%   The FUNCHEB class supports construction via interpolation at first- and
%   second-kind Chebyshev points with the classes FUNCHEB1 and FUNCHEB2
%   respectively. The default procedure is to use 2nd-kind points, but this can
%   be overwritted with the preferences PREF = FUNCHEB('tech', cheb1').
%
%   FUNCHEB.CONSTRUCTOR(VALUES, VSCALE, HSCALE, PREF) returns a FUNCHEB object
%   which interpolates the data in the columns of VALUES on a Chebyshev grid.
%   Whether this grid is of first- or second-kind points is determined by
%   PREF.FUNCHEB.tech, as above. FUNCHEB.CONSTRUCTOR({VALUES, COEFFS} allows for
%   the corresponding Chebyshev coefficients to be passed also.
%
% Examples:
%   % Basic construction:
%   f = funcheb.constructor(@(x) sin(x))
%
%   % Construction with preferences:
%   p = funcheb.pref('tech', 'cheb2'); % See help('funcheb.pref') for details.
%   f = funcheb.constructor(@(x) cos(x), p)
%
%   % Vector-valued construction:
%   f = funcheb.constructor(@(x) [sin(x), cos(x), exp(x)])
%
% See also FUNCHEB.PREF, HAPPINESSCHECK, FUNCHEB1, FUNCHEB2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB Class Description:
%
% The FUNCHEB class is an abstract class for representations of smooth functions
% on the interval [-1,1] via interpolated function values at Chebyshev points
% and coefficients of the corresponding first-kind Chebyshev series expansion.
%
% There are two main instances on the FUNCHEB class; FUNCHEB1 and FUNCHEB,
% which interpolate on Chebyshev grids of the first and second kind
% respectively. Note that although they use different Chebyshev grids in 'value'
% space, their coefficients are always from an expansion in first-kind Chebyshev
% polynomials (i.e., those usualy denoted by $T_k(x)$).
%
% The decision to use funcheb1 or funcheb2 is decided by the funcheb.pref.tech
% property, which should be either of the strings 'cheb1' or 'cheb2'.
%
% Class diagram: [<<smoothfun>>] <-- [<<FUNCHEB>>] <-- [funcheb1]
%                                                  <-- [funcheb2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of FUNCHEB objects.
    properties ( Access = public )
        
        % Values (stored in order from left to right) of FUNCHEB at Chebyshev
        % points. The particular Chebyshev points used depend on the instance of
        % the concrete class (1st kind for FUNCHEB1 and 2nd kind for FUNCHEB2).
        % For vectorised FUNCHEB objects, each column represents the
        % interpolated values of a single function.
        values % (nxm double)
        
        % Coefficients in 1st-kind Chebyshev series expansion of the FUNCHEB on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last. For vectorised FUNCHEB objects,
        % each column represents the coefficients of a single function.
        coeffs % (nxm double)
        
        % Vertical scale of the FUNCHEB. Typically the is the magnitude of the
        % largest value sampled from the given operator OP during the
        % construction process. It is updated via subsequent FUNCHEB operations
        % in a natural way.
        vscale = 0 % (1xm double >= 0)
        
        % Horizontal scale of the FUNCHEB. Although FUNCHEB objects have in
        % principle no notion of horizontal scale invariance (since they always
        % live on [-1,1]), the input OP has been implicitly mapped. HSCALE is
        % then used to enforce horizontal scale invariance in construction and
        % other subsequent operations that require it. It defaults to 1, and is
        % never updated.
        hscale = 1 % (scalar > 0)        
        
        % Boolean value designating whether the FUNCHEB is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)
        
        % Happiness level to which the FUNCHEB was constructed, or a rough
        % accuracy estimate of subsequent operations. See HAPPINESSCHECK.m for
        % full documentation.
        epslevel % (double >= 0)
    end
    
    %% CLASS CONSTRUCTOR:
    methods (Static)
        function obj = constructor(op, vscale, hscale, pref)
            % Constructor for the FUNCHEB class.
            
            % We can't return an empty FUNCHEB, so pass an empty OP down.
            if ( nargin == 0  )
                op = [];
            end
            
            % Define vscale if none given:
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            % Define vscale if none given:
            if ( nargin < 3 || isempty(hscale) )
                hscale = 1;
            end
            
            % Obtain preferences.
            if ( nargin == 2 && isstruct(vscale) )
                % vscale was actually a preference.
                pref = funcheb.pref(vscale);
                vscale = 0;
                hscale = 1;
            elseif ( nargin == 3 && isstruct(hscale) )
                % hscale was actually a preference.
                pref = funcheb.pref(hscale);
                hscale = 1;
            elseif ( nargin < 4 )
                % Create:
                pref = funcheb.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed.
                pref = funcheb.pref('eps', pref);
            else
                % Merge:
                pref = funcheb.pref(pref);
            end

            % Call the relevent constructor
            if ( strcmpi(pref.funcheb.tech, 'cheb1') )
                % Merge preferences:
                pref = funcheb1.pref(pref, pref.funcheb);
                % Construct:
                obj = funcheb1(op, vscale, hscale, pref);
            else
                % Merge preferences:
                pref = funcheb2.pref(pref, pref.funcheb);
                % Construct:
                obj = funcheb2(op, vscale, hscale, pref);
            end
            
        end
    end
    
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY THIS CLASS. 
    methods (Abstract)
        
        %
        
    end
    
    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods (Abstract, Static)
        
        % Compute Chebyshev barycentric weights.
        w = barywts(n)
        
        % Convert values to coefficients.
        coeffs = chebpoly(values)
        
        % Convert coefficients to values.
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points (x) and optionally quadrature (w) and
        % barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Extrapolate (for NaNs / endpoints).
        [values, maskNaN, maskInf] = extrapolate(values)
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Refinement function for FUNCHEB construction. (Evaluates OP on grid)
        [values, opints, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Make a FUNCHEB. (Constructor shortcut)
        f = make(varargin);
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods 
        
        %
        
    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        coeffs = alias(coeffs, m)
        
        % Evaluate a Chebyshev polynomial using barycentric interpolation.
        fx = bary(x, gvals, xk, vk, kind)

        % Evaluate a Chebyshev polynomial using Clenshaw's algorithm.
        out = clenshaw(x, coeffs)  
        
        % Test a sample evaluation of a FUNCHEB against and op evaluation.
        pass = sampleTest(op, f)
        
    end
    
end


