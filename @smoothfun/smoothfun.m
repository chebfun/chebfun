classdef smoothfun < onefun % (Abstract) 
%SMOOTHFUN   Approximate smooth functions on [-1,1]. 
%
%   Abstract (interface) class for approximating smooth functions on the
%   interval [-1,1].
%
% Constructor inputs:
%   SMOOTHFUN.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) constructs a SMOOTHFUN
%   object on the interval [-1,1] from the function handle OP. Currently the
%   only subclass of SMOOTHFUN is CHEBTECH, so SMOOTHFUN will call
%   CHEBTECH.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF2), where PREF2 is PREF merged
%   with the default CHEBTECH preferences.
%
% See also SMOOTHFUN.PREF, CHEBTECH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHFUN Class Description:
%
% The SMOOTHFUN class is an abstract class for representations of smooth
% functions on the interval [-1,1].
%
% Currently the only types of SMOOTHFUNs are CHEBTECH objects, which represent
% the smooth functions by Chebyshev interpolants.
%
% Class diagram: [<<onefun>>] <-- [<<SMOOTHFUN>>] <-- [<<chebtech>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Constructor for the SMOOTHFUN class.
    methods (Static)
        function obj = constructor(op, vscale, hscale, pref)
            
            % We can't return an empty SMOOTHFUN, so pass an empty OP down.
            if ( nargin == 0 )
                op = [];
            end
           
            % Define vscale if none given:
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            % Define hscale if none given:
            if ( nargin < 3 || isempty(hscale) )
                hscale = 1;
            end
            % Determine preferences if not given, merge if some are given:
            if ( nargin < 4 || isempty(pref) )
                pref = smoothfun.pref;
            else
                pref = smoothfun.pref(pref);
            end
            
            if ( strcmp(pref.smoothfun.tech, 'funqui') )
                op = funqui(op);
                pref.smoothfun.tech = smoothfun.pref('tech');
            end

            
            % Merge preferences:
            pref = chebtech.pref(pref, pref.smoothfun);
            % Call the CHEBTECH constructor
            obj = chebtech.constructor(op, vscale, hscale, pref);
            
        end
        
    end
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY SMOOTHFUN CLASS.
    methods ( Abstract = true )

    end

    %% ABSTRACT STATIC METHODS REQUIRED BY SMOOTHFUN CLASS.
    methods ( Abstract = true, Static = true )
        
    end
    
    %% Methods implimented by SMOOTHFUN class.
    methods 
        
    end
    
    %% Static methods implimented by SMOOTHFUN class.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Construct a rational interpolant to equispaced data.
        f = funqui(vals)
        
        % barycentric weights of Floater-Horman interpolant.
        w = fhBaryWts(n, d)
        
    end
    
end


function f = funqui(vals)
%FUNQUI    Rational interpolant of equispaced data.
%   F = FUNQUI(VALS) constructs a function handle F to a rational interpolant of
%   equispaced data in VALS in the interval [-1, 1]. It uses Floater-Hormann
%   interpolation with an adaptive choice of their blending degree d.

n = length(vals) - 1;

% Limit the maximal d to try depending on n:
if ( n < 60 )  
    maxd = min(n, 35);
elseif ( n < 100 )
    maxd = 30;
elseif ( n < 1000 )
    maxd = 25;
elseif ( n < 5000 )
    maxd = 20;
else
    maxd = 15;
end 

% Initialise:
errs = zeros(1, min(n, maxd) - 1);
x = linspace(-1, 1, n+1)';
xrm = x;
rmIndex = [2, n-1];
xrm(rmIndex) = [];
fvalsrm = vals;
fvalsrm(rmIndex) = [];

% Select a d:
if ( vals(rmIndex) < 2*eps*norm(vals, inf) ) 
    % This case fools funqui, so take a small d:
    dOpt = 4;
else
    % Find a near optimal d:
    for d = 0:min(n-2, maxd) 
        w = fhBaryWts(xrm, d);
        yyrm = chebtech.bary(x(rmIndex), fvalsrm, xrm, w);
        errs(d+1) = max( abs( yyrm - vals(rmIndex) ) );
        if ( errs(d+1) > 1000*min(errs(1:d+1)) )
            errs(d+2:end) = [];
            break
        end
    end
    % Find the index of the smallest error:
    [ignored, minInd] = min(errs); 
    dOpt = min(minInd) - 1;
end

% Compute FH weights:
w = fhBaryWts(x, dOpt);
% Create a function handle for the computed rational interpolant to the data:
f = @(zz) chebtech.bary(zz, vals, x, w);
end

function w = fhBaryWts(x, d) 
% Function for the computation of the FH weights.

n = length(x) - 1;
w = zeros(size(x));
for k = 1:n+1
   for m = k-d:k
      if ( m < 1 || m > n + 1 - d )
         continue
      end
      prod = 1;
      for j = m:m+d
         if ( j ~= k )
            prod = prod/( x(k) - x(j) );
         end
      end
      w(k) = w(k) + (-1)^m*prod;
   end
end
end
