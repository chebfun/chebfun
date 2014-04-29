classdef smoothfun < onefun % (Abstract) 
%SMOOTHFUN   Approximate smooth functions on [-1,1]. 
%   Abstract (interface) class for approximating smooth functions on the
%   interval [-1,1].
%
% Constructor inputs:
%   SMOOTHFUN.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) constructs a SMOOTHFUN
%   object on the interval [-1,1] from the function handle OP. Currently the
%   only subclass of SMOOTHFUN is CHEBTECH, so SMOOTHFUN will call
%   CHEBTECH.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF.TECHPREFS).
%
% See also CHEBTECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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
                pref = chebpref();
            else
                pref = chebpref(pref);
            end
            
            if ( strcmp(pref.tech, 'funqui') )
                op = funqui(op);
                if ( isfield(pref.techPrefs, 'funquiTech') )
                    pref.tech = pref.techPrefs.funquiTech;
                    pref.techPrefs = rmfield(pref.techPrefs, 'funquiTech');
                else
                    pref.tech = 'chebtech';
                end
            end

            % Call the CHEBTECH constructor
            obj = chebtech.constructor(op, vscale, hscale, pref.techPrefs);
            
        end
        
    end
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY SMOOTHFUN CLASS.
    methods ( Abstract = true )

    end

    %% ABSTRACT STATIC METHODS REQUIRED BY SMOOTHFUN CLASS.
    methods ( Abstract = true, Static = true )
        
    end
    
    %% Methods implemented by SMOOTHFUN class.
    methods 
        
    end
    
    %% Static methods implemented by SMOOTHFUN class.
    methods ( Static = true ) 
        
        % Construct a rational interpolant to equispaced data.
        f = funqui(vals)
        
        % barycentric weights of Floater-Horman interpolant.
        w = fhBaryWts(n, d)
        
    end
    
end

function f = funqui(vals)
%FUNQUI   Rational interpolant of equispaced data.
%   F = FUNQUI(VALS) constructs a function handle F to a rational interpolant of
%   equispaced data in VALS in the interval [-1, 1]. It uses Floater-Hormann
%   interpolation [Numer. Math. 107, 315-331 (2007)] with an adaptive choice 
%   of their blending degree d.

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

% Take arbitary linear combination of the columns for array-valued construction:
linComb = sin((1:size(vals, 2)).');
fvals = vals*linComb;
fvalsrm = fvals;
fvalsrm(rmIndex) = [];

% Select a d:
if ( norm(fvals(rmIndex), inf) < 2*eps*norm(fvals, inf) ) 
    % This case fools funqui, so take a small d:
    dOpt = 4;
else
    % Find a near optimal d:
    for d = 0:min(n - 2, maxd) 
        if ( d <= (n - 5)/2 )
            wl = abs(fhBaryWts(xrm, d, d + 2));
            wr = flipud(abs(fhBaryWts(flipud(xrm), d, d + 2)));
            w = [wl; wl(end)*ones(n - 5 - 2*d, 1); wr];
            w(1:2:end) = -w(1:2:end);
        else
            w = fhBaryWts(xrm, d);
        end
        yyrm = chebtech.bary(x(rmIndex), fvalsrm, xrm, w);
        errs(d+1) = max(abs(yyrm - fvals(rmIndex)));
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
if ( dOpt <= (n + 1)/2 ) 
    wl = abs(fhBaryWts(x, dOpt, dOpt + 1));
    w = [wl; wl(end)*ones(n - 1 - 2*dOpt, 1); flipud(wl)];
    w(1:2:end) = -w(1:2:end);
else
    w = fhBaryWts(x, dOpt);
end
% Create a function handle for the computed rational interpolant to the data:
f = @(zz) chebtech.bary(zz, vals, x, w);
end

function w = fhBaryWts(x, d, maxind) 
% Function for the computation of the FH weights.

n = length(x) - 1;

if ( nargin < 3 )
    maxind = n + 1;
end

w = zeros(min(n + 1, maxind), 1);

for k = 1:min(n + 1, maxind)
   for m = k-d:k
      if ( (m < 1) || (m > n + 1 - d) )
         continue
      end
      prod = 1;
      for j = m:m+d
         if ( j ~= k )
            prod = prod/(x(k) - x(j));
         end
      end
      w(k) = w(k) + (-1)^m*prod;
   end
end

end
