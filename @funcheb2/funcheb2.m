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
%   and the (column-wise) infinity norm of the sampled function values of OP,
%   and the fixed horizontal scale HSCALE. If not given, the VSCALE defaults to
%   0 initially, and HSCALE defaults to 1.
%
%   FUNCHEB2(OP, VSCALE, HSCALE, PREF) overrides the default behavior with that
%   given by the preference structure PREF. The constructor will also
%   accept inputs of the form FUNCHEB(OP, PREF), but this usage is not advised.
%   Similarly, one can pass FUNCHEB(OP, VSCALE, HSCALE, EPS), which is
%   equivalent to the call FUNCHEB(OP, VSCALE, HSCALE, FUNCHEB.PREF('eps',EPS)).
%
%   FUNCHEB2(VALUES, ...) returns a FUNCHEB2 object which interpolates the
%   values in the columns of VALUES at 2nd-kind Chebyshev points and
%   FUNCHEB2({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them from VALUES. If VALUES, is empty then the
%   FUNCHEB2 is constructed directly from the COEFFS.
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
% See also FUNCHEB, FUNCHEB2.pref, FUNCHEB2.chebpts, FUNCHEB2.happinesscheck, FUNCHEB2.refine.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB2 Class Description:
%
% The FUNCHEB2 class represents smooth functions on the interval [-1,1] using
% function values at 2nd-kind Chebyshev points and coefficients of the
% corresponding 1st-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 2nd-kind grid (see REFINE.m) until the
% representation is deemed 'happy' (see HAPPINESSCHECK.m). The resulting object
% can be used to evaluate and operate on the input function.
%
% More information can be found in the FUNCH class definition.
%
% Class diagram: [<<funcheb>>] <-- [FUNCHEB2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS IMPLEMENTED BY THIS M-FILE:
    methods
        
        function obj = funcheb2(op, vscale, hscale, pref)
            %Constructor for the FUNCHEB2 class.
            
            % Return an empty FUNCHEB2 on null input:
            if ( nargin == 0 || isempty(op) )
                return
            end
            
            % Define vscale if none given:
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            % Define hscale if none given:
            if ( nargin < 3 || isempty(hscale) )
                hscale = 1;
            end
            
            % Obtain preferences:
            if ( nargin == 2 && isstruct(vscale) )
                % vscale was actually a preference.
                pref = funcheb2.pref(vscale);
                vscale = 0;
                hscale = 1;
            elseif ( nargin == 3 && isstruct(hscale) )
                % hscale was actually a preference.
                pref = funcheb2.pref(hscale);
                hscale = 1;                
            elseif ( nargin < 4 )
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
            obj = populate(obj, op, vscale, hscale, pref);
            
            if ( obj.ishappy || isnumeric(op) || iscell(op) )
                % No need to error check if we are happy!
                return
            end
            
            % Check for NaNs: (if not happy)
            if ( pref.funcheb2.extrapolate )
                % Check for NaNs in interior only: (because extrapolate was on!)
                if ( any(any(isnan(obj.values(2:end-1,:)))) )
                    error('CHEBFUN:FUNCHEB2:constructor:naneval', ...
                        'Function returned NaN when evaluated.')
                end
                % We make sure not to return NaNs at +1 and -1.
                valuesTemp = funcheb2.extrapolate(obj.values);
                obj.values([1, end],:) = valuesTemp([1, end],:);
            elseif ( any(isnan(obj.values(:))) )
                % Here we throw an error if NaNs were encountered anywhere.
                error('CHEBFUN:FUNCHEB2:constructor:naneval2', ...
                    'Function returned NaN when evaluated.')
            end
            
        end
        
    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS:
    methods ( Static = true )
        
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
        
        % Extrapolate (for NaNs / endpoints).
        [values, maskNaN, maskInf] = extrapolate(values)
        
        % Make a FUNCHEB2. (Constructor shortcut)
        f = make(varargin);
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Refinement function for FUNCHEB2 construction. (Evaluates OP on grid)
        [values, points, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Test the FUNCHEB2 class.
        pass = test(varargin);
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS:
    methods
        
        % Compose two FUNCHEB2 objects or a FUNCHEB2 with a function_handle.
        h = compose(f, op, g, pref)
        
        % Get method.
        val = get(f, prop);
        
        % Set method.
%         f = set(f, prop, val); % [TODO]: Do we actually need a set method?
        
    end
    
end
