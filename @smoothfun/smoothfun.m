classdef smoothfun < onefun % (Abstract) 
%SMOOTHFUN   Approximate smooth functions on [-1,1]. 
%   Abstract (interface) class for approximating smooth functions on the
%   interval [-1,1].
%
% Constructor inputs:
%   SMOOTHFUN.CONSTRUCTOR(OP, DATA, PREF) constructs a SMOOTHFUN object on the
%   interval [-1,1] from the function handle OP using the data given in the
%   DATA structure and the preferences in PREF.
%
% See also ONEFUN, CHEBTECH, TRIGTECH.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHFUN Class Description:
%
% The SMOOTHFUN class is an abstract class for representations of smooth
% functions on the interval [-1,1].
%
% SMOOTHFUNs can be either TRIGTECH or CHEBTECH.
%
% Class diagram: [<<ONEFUN>>] <-- [<<SMOOTHFUN>>] <-- [TRIGTECH]
%                                                 <-- [<<CHEBTECH>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function obj = constructor(op, data, pref)
            
            % Parse inputs.
            if ( nargin == 0 )
                % We can't return an empty SMOOTHFUN, so pass an empty OP down.
                op = [];
            end

            if ( (nargin < 2) || isempty(data) )
                data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
            end

            % Deal with FUNQUI:
            if ( strcmp(pref.tech, 'funqui') )
                op = funqui(op);
                if ( isfield(pref.techPrefs, 'funquiTech') )
                    pref.tech = pref.techPrefs.funquiTech;
                    pref.techPrefs = rmfield(pref.techPrefs, 'funquiTech');
                else
                    pref.tech = @chebtech2;
                end
            end
            
            % Call the TECH constructor.
            obj = feval(pref.tech, op, data, pref.techPrefs);

        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false ) 
        
        function out = isPeriodicTech(f)
        %ISPERIODICTECH   Test if a SMOOTHFUN is constructed with a
        %basis of periodic functions. 
        %   Individual techs override this function as necessary.
            
            % Returns 0 by default.
            out = 0;
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Class-related functions: private utilities for this m-file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = funqui(vals)
%FUNQUI   Rational interpolant of equispaced data.
%   F = FUNQUI(VALS) constructs a function handle F to a rational interpolant of
%   equispaced data in VALS in the interval [-1, 1]. It uses Floater-Hormann
%   interpolation [Numer. Math. 107, 315-331 (2007)] with an adaptive choice 
%   of their blending degree d.

n = size(vals,1) - 1;

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
if ( n > 2 )
    rmIndex = [2, n-1];
else
    rmIndex = [];
end
xrm(rmIndex) = [];

% Take arbitrary linear combination of the columns for array-valued
% construction:
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
        yyrm = bary(x(rmIndex), fvalsrm, xrm, w);
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
    wm = repmat(wl(end), n - 1 - 2*dOpt, 1);
    w = [ wl ; wm ];
    w = w(1:ceil((n+1)/2));
    if ( mod(n, 2) )
        w = [ w ; flipud(w) ];
    else
        w = [ w(1:end-1) ; flipud(w) ];
    end
    w(1:2:end) = -w(1:2:end);
else
    w = fhBaryWts(x, dOpt);
end

% Create a function handle for the computed rational interpolant to the data:
f = @(zz) bary(zz, vals, x, w);

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
