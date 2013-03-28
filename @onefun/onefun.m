classdef onefun % (Abstract) 
%ONEFUN     Approximate smooth functions on [-1,1]. 
%   Abstract (interface) class for approximating functions on the interval
%   [-1,1].
%
% Constructor inputs:
%   ONEFUN.CONSTRUCTOR(OP) constructs a ONEFUN object from the function handle
%   OP. OP should be vectorised (i.e., accept a vector input) and ouput a vector
%   of the same length. Most ONEFUN objects allow for vectorised construction
%   (i.e., of multi-valued function), in which case OP should accept a vector of
%   length N and return a matrix of size NxM.
%
%   ONEFUN.CONSTRUCTOR(OP, VSCALE, HSCALE) constructs a ONEFUN with 'happiness'
%   relative to the maximum of the given vertical scale VSCALE (which is updated
%   by the infinity norm of the sampled function values of OP during
%   construction), and the fixed horizontal scale HSCALE. If not given, the
%   VSCALE defaults to 0 initially, and HSCALE defaults to 1.
%
%   ONEFUN.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) overrides the default behavior
%   with that given by the preference structure PREF. The constructor will also
%   accept inputs of the form ONEFUN(OP, PREF), but this usage is not advised.
%   Similarly, one can pass ONEFUN(OP, VSCALE, HSCALE, EPS), which is equivalent
%   to the call ONEFUN(OP, VSCALE, HSCALE, ONEFUN.PREF('eps',EPS))
%
%   ONEFUN.CONSTRUCTOR(VALUES, VSCALE, HSCALE, PREF) returns a ONEFUN object
%   which interpolates the values in the columns of VALUES. The points at which
%   this interpolation occurs is defined by PREF.ONEFUN.TECH.
%
% See also ONEFUN.pref, ONEFUN.chebpts, onefun.

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
        function obj = constructor(op, vscale, hscale, pref)
            
            % We can't return an empty ONEFUN, so pass an empty OP down.
            if ( nargin == 0 )
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
                pref = onefun.pref(vscale);
                vscale = 0;
                hscale = 1;
            elseif ( nargin == 3 && isstruct(hscale) )
                % hscale was actually a preference.
                pref = onefun.pref(hscale);
                hscale = 1;
            elseif ( nargin < 4 )
                % Create:
                pref = onefun.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed.
                pref = onefun.pref('eps', pref);
            else
                % Merge:
                pref = onefun.pref(pref);
            end

            % Call the relevent constructor:
            if ( isa(op, 'onefun') )
                % OP is already a ONEFUN!
                obj = op;
                
            elseif ( pref.onefun.blowup )
                % BLOWUP mode; Call SINGFUN.
                
                % Merge preferences:
                pref = singfun.pref(pref, pref.onefun);
                exponents = pref.singfun.exponents;
                
                % Call singfun constructor:
                obj = singfun(op, vscale, hscale, exponents, pref);
                
                % Return just a fun if no singularities found:
                if ( ~any(obj.exps) )
                    obj = obj.fun; 
                end 
                
            else
                % STANDARD mode; Call SMOOTHFUN.
                
                % Merge preferences:
                pref = smoothfun.pref(pref, pref.onefun);
                
                % Call SMOOTHFUN constructor:
                obj = smoothfun.constructor(op, vscale, hscale, pref);
                
            end
        
        end
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
    end
    
end
