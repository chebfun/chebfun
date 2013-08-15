function prefs = pref(varargin)
%PREF   Preference settings for CHEBFUN.
%   CHEBFUN.PREF(PREFNAME) returns the value corresponding to the preference
%   named in the string PREFNAME.
%
%   P = CHEBFUN.PREF returns a structure P with a field P.CHEBFUN which contains
%   the default CHEBFUN preferences as fields/values. This structure may be
%   passed to the CHEBFUN constructor.
%
%   P = CHEBFUN.PREF(P) will check to see whether the input preference structure
%   P already has a CHEBFUN field. If it does not, one is appended.
%
%   P = CHEBFUN.CHEBFUN('PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) returns the
%   same structure as above, but with the default CHEBFUN preferences
%   'PREFNAME1', 'PREFNAME2', replaced by the values in VAL1, VAL2, etc.
%
%   P = CHEBFUN.PREF(P, 'PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) appends a
%   CHEBFUN preference field to P if required, and modifies the FUN properties
%   'PREFNAME1' and 'PREFNAME2'.
%
%   Note that no checking of either the input PREFNAMEs or VALs takes place.
%
%   CHEBFUN PREFERENCES (case sensitive)
%
%    domain
%      [-1,1]      - Default domain for FUN construction.
%
%     eps          - Relative tolerance used in construction and subsequent
%      [2^-52]       operations.
%
%    extrapolate   - Extrapolation at endpoints.
%     [0, 'off']     If 'on', function values at endpoints maybe extrapolated
%      1, 'on'       from interior values rather than sampled. Extrapolation is
%                    used independently of this option if a function evaluation
%                    returns NaN. In some cases, however, function values at
%                    end points may be inaccurate or undefined, and enabling
%                    extrapolation may be helpful.
%
%    minsamples    - Minimum number of points used by the constructor. The 
%       [9]          constructed CHEBFUN might be shorter as a result of
%                    the SIMPLIFY() command. Must be of the form 2^n+1.
%
%    maxdegree     - Maximum degree used by constructor in SPLITTING 'OFF' mode.
%      [2^16]
%
%    splitting     - Domain splitting option.
%     [0, 'off']     If 'on', breakpoints between funs may be introduced where a
%      1, 'on        discontinuity in a function or a low-order derivative is
%                    detected, or if a global polynomial representation will be
%                    too long. If 'off', breakpoints will be introduced only at
%                    points where discontinuities are being created, e.g., by
%                    ABS(F) at points where a CHEBFUN F passes through zero.
%
%    splitdegree   - Maximum degree used by constructor in SPLITTING 'ON' mode.
%       [128] 
%
%     maxlength    - Maximum total number of points used in SPLITTING 'ON' mode.
%       [6000]
%
% See also CHEBFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% blowup - Blowup option:
%        BLOWUP=0: bounded functions only
%        BLOWUP=1: poles are permitted (integer order)
%        BLOWUP=2: blowup of integer or non-integer orders permitted (experimental)

classname = 'chebfun';

% Has a preference structure been passed?
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    prefs = varargin{1};
    varargin(1) = [];
else                                                % No
    prefs = struct();
end

% If it was, did it have a CHEBFUN field?
if ( isfield(prefs, classname) )  % It does, so either:
    if ( numel(varargin) == 0 )
        return                    % a) No prefs to change, return.
    else
        p = prefs.(classname);    % b) Grab chebfun prefs.
    end
else                              % No chebfun prefs found, so make some:
    %p.blowup      = false;
    p.eps         = 2^-52;
    p.extrapolate = true;
    p.domain      = [-1, 1];
    p.splitting   = false;
    % TODO: Chebfun should not know about 'degree'.
    p.maxdegree   = 65536;
    p.maxlength   = 6000;
    p.splitdegree = 128;
    p.tech        = 'chebtech2';
end
% p is now the preference substructure relating to the CHEBFUN class.

if ( isfield(prefs,'misc') )
    q = prefs.misc;
else
    q = struct();
end
% q is now the preference substructure for MISC preferences

for names = fieldnames(q)'
    if isfield(p,names)
        field = names{:};
        p.(field) = q.(field);
        q = rmfield(q,field);
    end
end

% Two preference structures were passed. Copy matching fieldnames.
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

% A single property has been queried, so return this property
if ( numel(varargin) == 1 )
    prefs = p.(varargin{1});
    return
end

% Property names have been passed, so alter/add CHEBFUN/MISC properties.
for k = 1:2:numel(varargin)
    if ( strcmp(varargin{2},'off') )
        varargin{2} = false;
    elseif ( strcmp(varargin{2},'on') )
        varargin{2} = true;
    end
    if ( isfield(p, varargin{k}) )
        p.(varargin{k}) = varargin{k+1};
    else
        q.(varargin{k}) = varargin{k+1};
    end
end

% Append CHEBFUN preferences to the preference structure prefs for output.
prefs.(classname) = p;

% Append MISC preferences to the preference structure prefs for output.
prefs.misc = q;

