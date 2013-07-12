classdef singfun
    
%SINGFUN Class for functions with singular endpoint behavior.
%   TODO: User documentation

%% SINGFUN class description
%
% The singfun class represents a function of the form 
%
% \[ f(x) = s(x) (1+x)^\alpha (1-x)^\beta \]
%
% on the interval $[-1,1]$. The exponents $\alpha$ and $\beta$ are assumed
% to be real and negative. The constructor is supplied with a handle that evaluates the
% function $f$ at any given points. However, endpoint values will not be
% sampled, due to the likelihood of Inf and NaN results.
%
% Ideally, the "smooth" function $s$ is analytic, or at least much more
% compactly represented than $f$ is. The resulting object can be used to
% evaluate and operate on the function $f$. If $\alpha$ and $\beta$ are
% unknown at the time of construction, the constructor will try to
% determine appropriate (nonpositive) values automatically by sampling
% the function handle. Note, however, that this process is not 
% completely robust, and
% the singularity terms in general do not perfectly factor out singular
% behavior. The constructor can be forced to consider only integer
% exponents.
%
% Multiplication and division are as good as the corresponding operations
% on the smooth part. Addition and subtraction are much less reliable, as
% the sum of two singfuns with different exponents is not necessarily a
% singfun, nor a smooth function. If all but integer exponents can be
% factored out of the summands, the process is fine, but in other
% circumstances the process may throw an error.
 
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of SINFGUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        smoothPart  % (smoothfun)
        
        % Exponents of the singularities at the two endpoints.
        exponents   % (1x2 double)
        
        % A cell array telling the type of singularity at the endpoints.
        singType   % (1x2 cell)
        
        % A logical array indicating which ends are singular.
        isSingEnd   % (1x2 logical)        
    end
    
    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = singfun(op, singFlag, type, pref)           
            %%
            % no input arguments: return an empty object               
            if ( nargin == 0 )   
                obj.smoothPart = [];
                obj.exponents = [];
                obj.singType = {};
                obj.isSingEnd = [];
                return
            end
            
            if ( nargin == 1 )
                % only operator passed, assume a pole at each end point
                obj.singType = {'pole', 'pole'};
                obj.isSingEnd = [1, 1];
            else
                % copy the information given about singularities in the current object
                obj.singType = type;
                obj.isSingEnd = singFlag;
            end
                                   
            % Determine preferences if not given, merge if some are given:
            if ( nargin < 4 || isempty(pref) )
                pref = singfun.pref;
            else        
                pref = singfun.pref(pref);
            end
            
            %%
            % Determine and factor out singular terms.
            obj.exponents = singfun.findSingExponents(op, obj.isSingEnd, obj.singType, pref);
            
            % update ISSINGEND and SINGTYPE based on EXPONENTS
            tol = pref.singfun.eps;
            if ( abs(obj.exponents(1)) < 100*tol )
                % if the singularity exponent is below the tolerance level
                % remove the singularity
                obj.isSingEnd(1) = 0;
                obj.singType{1} = 'none';
            end            
            if ( abs(obj.exponents(2)) < 100*tol )
                % if the singularity exponent is below the tolerance level
                % remove the singularity
                obj.isSingEnd(2) = 0;
                obj.singType{2} = 'none';
            end
                                  
            % We check for three cases to avoid doubly nesting functions
            % when both exponents are present. 
            if ( all(abs(obj.exponents) > 100*tol ) )         
                % left and right terms
                op = @(x) op(x)./((1+x).^(obj.exponents(1)).*(1-x).^(obj.exponents(2)));
            elseif ( abs(obj.exponents(1)) > 100*tol )   
                % left only
                op = @(x) op(x)./(1+x).^(obj.exponents(1));
            elseif ( abs(obj.exponents(2)) > 100*tol )    
                % right only
                op = @(x) op(x)./(1-x).^(obj.exponents(2));
            end
            
            % Construct the smooth part of the SINGFUN object.
            % [TODO]: This will be replaced by the SMOOTHFUN constructor
            prefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
            vscale = [];
            hscale = [];
            obj.smoothPart = chebtech.constructor(op, vscale, hscale, prefs);
        end
    end
    
    %% Public methods defined in external files.
    methods ( Access = public )
        y = feval(f,x)
    end
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods ( Static = true )
        
        exponents = findSingExponents( op, isSingEnd, singType, pref )
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)

    end

end
    