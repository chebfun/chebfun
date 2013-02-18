classdef fun2 < smoothfun
%FUN2   Approximate smooth functions on [-1,1] with Chebyshev interpolants. 
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values at 2nd-kind Chebyshev points and coefficients of the
%   corresponding first-kind Chebyshev series expansion.
%
% Constructor inputs:
%   FUN2(OP) constructs a FUN2 object from the function handle OP. OP should be
%   vectorized (i.e., accept a vector input), and ouput a vector of the same
%   length.
%
%   FUN2(OP, VSCALE) constructs a FUN2 with happiness relative to the given
%   vertical (VSCALE) scale. If not given, thise value defaults 0.
%
%   FUN2(OP, VSCALE, PREF) overrides the default behavior with that given by the
%   chebpref object PREF.
%
% Examples:
%   % Basic construction
%   f = fun2(@(x) sin(x))
%
%   % Construction with preferences
%   prefs = fun2.pref('sampletest', 0); % See help('fun2.pref') for details
%   f = fun2(@(x) sin(x), [], [], prefs)
%
% See also FUN2.pref, FUN2.chebpts, smoothfun.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN2 Class Description:
%
% The FUN2 class represents smooth functions on the interval $[-1,1]$ using
% function values at 2nd-kind Chebyshev points and coefficients of the
% corresponding first-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 2nd-kind grid (see refinementStrategy.m ?)
% until the representation is deemed 'happy' (see below, as well as ishappy.m
% and simplify.m). The resulting object can be used to evaluate and operate on
% the input function.
%
% [TODO] Explain VSCALE, HLEVEL, and how the happiness check works
% (sampletest.m).
%
% [TODO] Explain extrapolation.
%
% [TODO] Experiment with vectorised functions.
%
% Class diagram: [smoothfun] <-- [fun2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Properties of FUN2 objects.
    properties ( Access = public )
        
        % Values of FUN2 at 2nd-kind Chebyshev points. Values are stored in
        % order from left to right (i.e., values(1) = feval(FUN2,-1)).
        values % (nx1 double)
        
        % Coefficients in first-kind Chebyshev series expansion of the FUN2 on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last.
        coeffs % (nx1 double)
        
        % Vertical scale of the FUN2. [TODO]: more detail.
        vscale % (double >= 0)
        
        % Happiness level to which the FUN2 was constructed. [TODO]: more
        % detail.
        hlevel % (double >= 0)
    end
    
    %% Methods implemented by this m-file.
    methods
        
        function obj = fun2(op, vscale, pref)
            % Constructor for the FUN2 class.
            
            % Return an empty FUN2 on null input.
            if ( nargin == 0 || isempty(op) )
                return
            end
            
            % Define vscale if none given.
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            
            % Obtain preferences if none given.
            if ( nargin < 3 )
                % Create
                pref = fun2.pref;
            else
                % Merge
                pref = fun2.pref(pref);
            end
            
            % Force nonadaptive construction.
            if ( isfield(pref.fun2, 'n') && ~isnan(pref.fun2.n) )
                % Note that pref.fun2.n = NaN will be adaptive
                op = feval(op, fun2.chebpts(pref.n));
            end
            
            % Actual construction takes place here.
            if ( ~isnumeric(op) )
                % Adaptive contruction.
                [values, coeffs, vscale, hlevel] = ...
                    fun2.constructor(op, vscale, pref); %#ok<*PROP>
            else
                % Nonadaptive contruction.
                values = op;
                coeffs = fun2.chebpoly(values);
                vscale = norm(values, inf);
                hlevel = -1; % [TODO] What should this return?
            end
            
            % Assign to FUN2 object.
            obj.values = values;
            obj.coeffs = coeffs;
            obj.vscale = vscale;
            obj.hlevel = hlevel;
            
        end
        
    end
    
    %% Static methods implemented by FUN2 class.
    methods ( Static = true )
        
        % The main constructor
        [values, coeffs, vscale, hlevel] = ...
            constructor(op, vscale, pref)
        
        % Happiness check
        [hlevel, values, coeffs] = ...
            ishappy(op, values, coeffs, vscale, pref)
        
        % Simplify Chebyshev coefficients
        [values, coeffs, hlevel] = ...
            simplify(values, coeffs, vscale, pref)
        
        % Extrapolate (for NaNs / endpoints)
        [values, maskNaN, maskInf] = extrapolate(values, vscale, pref)
        
        % Convert values to coefficients
        coeffs = chebpoly(values)
        
        % Convert coefficients to values
        values = chebpolyval(coeffs)
        
        % Compute chebyshev points of first kind (x) and optionally
        % quadrature (w) and barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Evaluate a Chebyshev polynomial using Clenshaw's algorithm.
        out = clenshaw(coeffs, x)
        
        % Alias Chebyshev coefficients
        coeffs = alias(coeffs, m);
        
        % Test the FUN2 class.
        pass = test(varargin);
        
    end
    
    %% Private static methods implemented by FUN2 class.
    methods ( Access = private, Static = true )
        
        % Define the refinement strategy used during construction. This can be
        % overridden as a FUN2 pref.
        nn = refinementStrategy(pref)
        
    end
    
end
