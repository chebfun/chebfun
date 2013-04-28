classdef chebtech1 < chebtech
%CHEBTECH1   Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Class for approximating smooth functions on the interval [-1,1]
%   using function values at 1st-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   CHEBTECH1(OP) constructs a CHEBTECH1 object from the function handle OP. OP
%   should be vectorized (i.e., accept a vector input) and ouput a vector of
%   the same length. CHEBTECH1 objects allow for vectorised construction (i.e.,
%   of array-valued function), in which case OP should accept a vector of length
%   N and return a matrix of size NxM, where M is number of columns of the multi
%   -valued function.
%
%   CHEBTECH1(OP, VSCALE) constructs a CHEBTECH1 with 'happiness' (see
%   CHEBTECH.HAPPINESSCHECK) relative to the maximum of the given vertical scale 
%   (VSCALE) and the (column-wise) infinity norm of the sampled function values 
%   of OP. If not given, the VSCALE defaults to 0 initially.
%
%   CHEBTECH1(OP, VSCALE, HSCALE) uses a 'happiness' to both the vertical scale
%   VSCALE (as above) and the horizontal scale HSCALE. If not given (or given as
%   empty), this defaults to 1.
%
%   CHEBTECH1(OP, VSCALE, HSCALE, PREF) overrides the default behavior with that
%   given by the preference structure PREF.
%
%   CHEBTECH1(VALUES, ...) returns a CHEBTECH1 object which interpolates the
%   values in the columns of VALUES at 1st-kind Chebyshev points and
%   CHEBTECH1({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them. If COEFFS are passed, the resulting
%   CHEBTECH1 is always deemed 'happy'.
%
% Examples:
%   % Basic construction:
%   f = chebtech1(@(x) sin(x))
%
%   % Construction with preferences:
%   p = chebtech.pref('sampletest', 0); % See help('chebtech.pref') for details
%   f = chebtech1(@(x) sin(x), [], [], p)
%
%   % Vector-valued construction:
%   f = chebtech1(@(x) [sin(x), cos(x), exp(x)])
%
% See also CHEBTECH, CHEBTECH.PREF, CHEBPTS, HAPPINESSCHECK, REFINE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBTECH1 Class Description:
%
% The CHEBTECH1 class represents smooth functions on the interval [-1,1] using
% function values at 1st-kind Chebyshev points and coefficients of the
% corresponding 1st-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 1st-kind grid (see REFINE.m) until the
% representation is deemed 'happy' (see HAPPINESSCHECK). The resulting
% object can be used to evaluate and operate on the input function.
%
% More information can be found in the CHEBTECH class definition.
%
% Class diagram: [<<CHEBTECH>>] <-- [CHEBTECH1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS IMPLEMENTED BY THIS M-FILE:
    methods
        
        function obj = chebtech1(op, vscale, hscale, pref)
            % Constructor for the CHEBTECH1 class.
            
            % Return an empty CHEBTECH1 on null input:
            if ( (nargin == 0) || isempty(op) )
                return
            end
           
            % Define vscale if none given:
            if ( (nargin < 2) || isempty(vscale) )
                vscale = 0;
            end
            % Define hscale if none given:
            if ( (nargin < 3) || isempty(hscale) )
                hscale = 1;
            end
            % Determine preferences if not given, merge if some are given:
            if ( (nargin < 4) || isempty(pref) )
                pref = chebtech.pref;
            else
                pref = chebtech.pref(pref);
            end

            % Force nonadaptive construction if PREF.CHEBTECH.N is numeric:
            if ( ~isempty(pref.chebtech.n) && ~isnan(pref.chebtech.n) )
                % Evaluate op on the Chebyshev grid of given size:
                op = feval(op, chebtech1.chebpts(pref.chebtech.n));
            end
            
            % Actual construction takes place here:
            obj = populate(obj, op, vscale, hscale, pref);
            
            if ( obj.ishappy || isnumeric(op) || iscell(op) )
                % No need to error check if we are happy:
                return
            end
            
            % Check for NaNs (if not happy):
            if ( any(isnan(obj.values(:))) )
                % Here we throw an error if NaNs were encountered anywhere.
                error('CHEBFUN:CHEBTECH1:constructor:naneval', ...
                    'Function returned NaN when evaluated.')
            end
            
        end
        
    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS:
    methods ( Static = true )
        
        % Aliasing:
        coeffs = alias(coeffs, m)
        
        % Evaluate a Chebyshev interpolant using proper barycentric formula:
        out = bary(x, values)
        
        % Compute Chebyshev barycentric weights:
        w = barywts(n)

        % Convert values to coefficients:
        coeffs = chebpoly(values)
        
        % Convert coefficients to values:
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points (x) and optionally quadrature (w)
        % and barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Make a CHEBTECH1 (constructor shortcut):
        f = make(varargin);

        % Compute Chebyshev quadrature weights:
        w = quadwts(n)
        
        % Refinement function for CHEBTECH1 construction (evaluates OP on grid):
        [values, points, giveUp] = refine(op, values, pref)
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS:
    methods
        
        % Compose two CHEBTECH1 objects or a CHEBTECH1 with a function handle:
        h = compose(f, op, g, pref)
        
        % Get method:
        val = get(f, prop);

    end
    
end
