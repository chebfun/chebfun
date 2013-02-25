classdef (Abstract) onefun
%ONEFUN     Approximate smooth functions on [-1,1]. 
%   Abstract (interface) class for approximating functions on the interval
%   [-1,1].
%
% Constructor inputs:
%   ONEFUN(OP) constructs a ONEFUN object from the function handle OP. OP
%   should be vectorised (i.e., accept a vector input) and ouput a vector of the
%   same length. Most ONEFUN objects allow for vectorised construction (i.e.,
%   of multi-valued function), in which case OP should accept a vector of length
%   N and return a matrix of size NxM.
%
%   ONEFUN.CONSTRUCTOR(OP, VSCALE) constructs a ONEFUN with 'happiness'
%   relative to the maximum of the given vertical scale VSCALE and the infinity
%   norm of the sampled function values of OP. If not given, the VSCALE defaults
%   to 0 initially.
%
%   ONEFUN.CONSTRUCTOR(OP, VSCALE, PREF) overrides the default behavior with
%   that given by the preference structure PREF. The constructor will also
%   accept inputs of the form ONEFUN(OP, PREF), but this usage is not
%   advised. Similarly, one can pass ONEFUN(OP, VSCALE, EPS), which is
%   equivalent to the call ONEFUN(OP, VSCALE, ONEFUN.PREF('eps',EPS))
%
%   ONEFUN.CONSTRUCTOR(VALUES, VSCALE, PREF) returns a ONEFUN object which
%   interpolates the values in the columns of VALUES. The points at which this
%   interpolation occurs is defined by PREF.ONEFUN.TECH.
%
% See also ONEFUN.pref, ONEFUN.chebpts, onefun.

% [TODO]: Document this file.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONEFUN Class Description:
%
% The ONEFUN class is an abstract class for representations of functions on the
% interval [-1,1].
%
% The current instances of ONEFUNs are SMOOTHFUNs and SINGFUNs. The former are
% used to represenet smooth functions on [-1,1}, whereas the latter are able to
% represent some forms of endpoint singularites. 
%
% Class diagram: [<<fun>>] <>-- [<<ONEFUN>>] <-- [<<smoothfun>>]
%                                            <-- [singfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    methods (Static)
        function y = constructor(op, vscale, hscale, pref)
            
            if ( nargin < 3 )
                pref = onefun.pref;
            else
                pref = onefun.pref(pref);
            end
            
            if ( nargin < 2 )
                vscale = 0;
            end

            % Call the relevent constructor
            if ( pref.onefun.blowup )
                
                pref = singfun.pref(pref, pref.onefun);
                
                % Call singfun constructor
                y = singfun(op, exponents, pref);
                
                % Return just a fun if no singularities found
                if ( ~any(y.exps) )
                    y = y.fun; 
                end 
                
            else
                
                pref = smoothfun.pref(pref, pref.onefun);
                
                % Call SMOOTHFUN constructor
                y = smoothfun.constructor(op, vscale, hscale, pref);
                
            end
        
        end
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
    end
    
end
