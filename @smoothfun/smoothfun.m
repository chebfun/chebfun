classdef (Abstract) smoothfun
%SMOOTHFUN   Approximate smooth functions on [-1,1]. 
%
%   Abstract (interface) class for approximating smooth functions on the
%   interval [-1,1].
%
% Constructor inputs:
%   SMOOTHFUN(OP) constructs a SMOOTHFUN object from the function handle OP. OP
%   should be vectorised (i.e., accept a vector input) and ouput a vector of the
%   same length. Most SMOOTHFUN objects allow for vectorised construction (i.e.,
%   of multi-valued function), in which case OP should accept a vector of length
%   N and return a matrix of size NxM.
%
%   SMOOTHFUN.CONSTRUCTOR(OP, VSCALE) constructs a SMOOTHFUN with 'happiness'
%   relative to the maximum of the given vertical scale VSCALE and the infinity
%   norm of the sampled function values of OP. If not given, the VSCALE defaults
%   to 0 initially.
%
%   SMOOTHFUN.CONSTRUCTOR(OP, VSCALE, PREF) overrides the default behavior with
%   that given by the preference structure PREF. The constructor will also
%   accept inputs of the form SMOOTHFUN(OP, PREF), but this usage is not
%   advised. Similarly, one can pass SMOOTHFUN(OP, VSCALE, EPS), which is
%   equivalent to the call SMOOTHFUN(OP, VSCALE, SMOOTHFUN.PREF('eps',EPS))
%
%   SMOOTHFUN.CONSTRUCTOR(VALUES, VSCALE, PREF) returns a SMOOTHFUN object which
%   interpolates the values in the columns of VALUES. The points at which this
%   interpolation occurs is defined by PREF.SMOOTHFUN.TECH.
%
% See also SMOOTHFUN.pref, SMOOTHFUN.chebpts, onefun.

% [TODO]: Document this file.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHFUN Class Description:
%
% The SMOOTHFUN class is an abstract class for representations of smooth
% functions on the interval [-1,1].
%
% Currently the only types of SMOOTHFUNs are funcheb obkects, which represent
% the smooth functions by Chebyshev interpolants.
%
% Class diagram: [<<onefun>>] <-- [<<SMOOTHFUN>>] <-- [<<funcheb>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.
    
    methods (Static)
        
        function obj = constructor(op, vscale, pref)
            % Constructor for the SMOOTHFUN class.
            
            % We can't return an empty SMOOTHFUN, so pass an empty OP down.
            if ( nargin == 0 )
                op = [];
            end
            
            % Define vscale if none given.
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            
            % Obtain preferences.
            if ( nargin == 2 && isstruct(vscale) )
                % vscale was actually a preference.
                pref = smoothfun.pref(vscale);
                vscale = 0;
            elseif ( nargin < 3 )
                % Create
                pref = smoothfun.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed
                pref = smoothfun.pref('eps', pref);
            else
                % Merge
                pref = smoothfun.pref(pref);
            end

            % Call the relevent constructor
            pref = funcheb.pref(pref, pref.smoothfun);
            obj = funcheb.constructor(op, vscale, pref);
            
        end
        
    end
    
    % Methods implimented by SMOOTHFUN class.
    methods 
        
    end
    
    % Static methods implimented by SMOOTHFUN class.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
    end
    
end


