function prefs = pref(varargin)
%PREF   Preference settings for CHEBTECH.
%   CHEBTECH.PREF(PREFNAME) returns the value corresponding to the preference
%   named in the string PREFNAME.
%
%   P = CHEBTECH.PREF returns a structure P with a field P.CHEBTECH which contains
%   the default CHEBTECH preferences as fields/values. This structure may be
%   passed to the CHEBTECH constructor.
%
%   P = CHEBTECH.PREF(P) will check to see whether the input preference structure
%   P already has a CHEBTECH field. If it does not, one is appended.
%
%   P = CHEBTECH.PREF('PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) returns the same
%   structure as above, but with the default CHEBTECH preferences 'PREFNAME1',
%   'PREFNAME2', replaced by the values in VAL1, VAL2, etc.
%
%   P = CHEBTECH.PREF(P, 'PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) appends a
%   CHEBTECH preference field to P if required, and modifies the CHEBTECH
%   properties 'PREFNAME1' and 'PREFNAME2'.
%
%   Note that no checking of either the input PREFNAMEs or VALs takes place.
%
%   CHEBTECH PREFERENCES (case sensitive)
%
%     tech         -  Select the type of Chebyshev grid on which the function
%                     is sampled.
%   ['chebtech2']  -  Use a grid of second-kind Chebyshev points.
%    'chebtech1'   -  Use a grid of first-kind Chebyshev points.
%
%     eps          -  Relative tolerance used in construction and subsequent
%      [2^-52]        operations. See CHEBTECH.HAPPINESSCHECK for more details.
%
%     extrapolate
%        true      -  Function values at endpoints may be extrapolated from
%                     interior values rather than sampled.
%       [false]    -  Do not extrapolate values at endpoints.
%
%     hscale       -  Horizontal scale. This preference can be used to ensure
%        [1]          horizontal scale invariance when using the CHEBTECH
%                     constructor to implicitly represent functions defined on
%                     domains other than [-1, 1].
%
%     minSamples   -  Minimum number of points used by the constructor. Should
%        [9]          be of the form 2^n+1 (if not, it is rounded as such).
%
%     maxSamples   -  Maximum number of points used by the constructor.
%      [2^16+1]   
%
%         n        -  Fixed number of points used by constructor. NaN allows
%       [NaN]         adaptive construction.
%
%     sampletest
%       [true]     -  Tests the function at one more arbitrary point to
%                     minimize the risk of missing signals between grid
%                     points.
%        false     -  Do not test.
%
%   refinementFunction - Define function for refining sample values.
%     ['nested']       - Use the default process (nested evaluation).
%      'resampling'    - Every function value is computed afresh as the
%                        constructor tries grids of size 9, 17, 33, etc.
%   function_handle    - A user-defined refinement. See REFINE.m
%
%   happinessCheck     - Define function for testing happiness.
%     ['classic']      - Use the default process from Chebfun V4.
%      'strict'        - Strict tolerance for coefficients.
%      'loose'         - A looser tolerance for coefficients.
%     function_handle  - A user defined happiness. See HAPPINESSCHECK.m
%
% See also CHEBTECH, CHEBTECH1, CHEBTECH2

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

classname = 'chebtech';

% Has a preference structure been passed?
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    prefs = varargin{1};
    varargin(1) = [];
else                                                % No
    prefs = struct();
end

% If it was, did it have a field for this class?
if ( isfield(prefs, classname) )  % It does, so either:
    if ( numel(varargin) == 0 )
        return                    % a) No props to change, return
    else
        p = prefs.(classname);    % b) Grab prefs for this class.
    end
else                              % No prefs found for this class, so make some:
    p.tech        = 'chebtech2';
    p.eps         = 2^-52;
    p.extrapolate = false;
    p.hscale      = 1;
    p.minSamples  = 9;    
    p.maxSamples  = 2^16+1;
    p.n           = NaN;
    p.sampletest  = true;
    p.refinementFunction = 'nested';
    p.happinessCheck = 'classic';
end
% p is now the preference substructure relating to the current class.

if ( isfield(prefs,'misc') ) 
    q = prefs.misc;
else
    q = struct();
end
% q is now the preference substructure for MISC preferences.

for names = fieldnames(q)'
    if isfield(p,names)
        field = names{:};
        p.(field) = q.(field);
        q = rmfield(q,field);
    end
end

% Two preference structures were passed. Copy matching fieldnames:
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    r = varargin{1};
    varargin(1) = [];
    for names = fieldnames(r).'
        if isfield(p,names)
            field = names{:};
            p.(field) = r.(field);
        end
    end
end

% A single property has been queried, so return this property:
if ( numel(varargin) == 1 )
    prefs = p.(varargin{1});
    return
end

% Property names have been passed, so alter/add <current class>/MISC properties.
for k = 1:2:numel(varargin)
    if ( isfield(p, varargin{k}) )
        p.(varargin{k}) = varargin{k+1};
    else
        q.(varargin{k}) = varargin{k+1};
    end
end

% Append current class preferences to the preference structure prefs for output:
prefs.(classname) = p;

% Append MISC preferences to the preference structure prefs for output:
prefs.misc = q;

end
