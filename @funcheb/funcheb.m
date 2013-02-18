classdef (Abstract) funcheb < smoothfun
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
%   of multi-valued function), in which case OP should accept a column vector of
%   length N and return a matrix of size NxM.
%
%   FUNCHEB.CONSTRUCTOR(OP, VSCALE) constructs a FUNCHEB with 'happiness'
%   relative to the maximum of the given vertical scale VSCALE and the infinity
%   norm of the sampled function values of OP. If not given, the VSCALE defaults
%   to 0 initially.
%
%   FUNCHEB.CONSTRUCTOR(OP, VSCALE, PREF) overrides the default behavior with
%   that given by the preference structure PREF. The constructor will also
%   accept inputs of the form FUNCHEB(OP, PREF), but this usage is not advised.
%   Similarly, one can pass FUNCHEB(OP, VSCALE, EPS), which is equivalent to the
%   call FUNCHEB(OP, VSCALE, FUNCHEB.PREF('eps',EPS))
%
%   FUNCHEB.CONSTRUCTOR(VALUES, VSCALE, PREF) returns a FUNCHEB object which
%   interpolates the values in the columns of VALUES.
%
% Examples:
%
% See also 

% [TODO]: Document this file.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB Class Description:
%
% The FUNCHEB class is an abstract class for representations of smooth functions
% on the interval [-1,1] via interpolated function values at Chebyshev points
% and coefficients of the corresponding first-kind Chebyshev series expansion.
%
% There are two main instances on the FUNCHEB class; FUNCHEB1 and FUNCHEB2,
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
    
    methods (Static)
        
        function obj = constructor(op, vscale, pref)
            % Constructor for the SMOOTHFUN class.
            
            % We can't return an empty FUNCHEB, so pass an empty OP down.
            if ( nargin == 0  )
                op = [];
            end
            
            % Define vscale if none given.
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            
            % Obtain preferences.
            if ( nargin == 2 && isstruct(vscale) )
                % vscale was actually a preference.
                pref = funcheb.pref(vscale);
                vscale = 0;
            elseif ( nargin < 3 )
                % Create
                pref = funcheb.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed
                pref = funcheb.pref('eps', pref);
            else
                % Merge
                pref = funcheb.pref(pref);
            end

            % Call the relevent constructor
            if ( strcmpi(pref.funcheb.tech, 'cheb1') )
                pref = funcheb1.pref(pref, pref.smoothfun);
                obj = funcheb1(op, vscale, pref);
            else
                pref = funcheb2.pref(pref, pref.smoothfun);
                obj = funcheb2(op, vscale, pref);
            end
            
        end
        
    end
    
    methods % Methods implimented by SMOOTHFUN class.
        
    end
    
    methods ( Static = true ) % Static methods implimented by SMOOTHFUN class.
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
    end
    
end


