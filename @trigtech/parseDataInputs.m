function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.
%   DATA = PARSEDATAINPUTS(DATA, PREF) parses inputs and options passed to the
%   constructor via its DATA argument and returns a new DATA structure with the
%   results, including defaults for any missing values.  The preferences in
%   PREF are used as necessary to determine the results.
%
%   The DATA argument to the constructor is a structure that holds all
%   additional parameters required for construction aside from the OP (the
%   function handle that is sampled or numeric data that is used to build the
%   TRIGTECH) and the constructor preferences.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isfield(data, 'vscale') || isempty(data.vscale) )
    data.vscale = 0;
end

if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    data.hscale = 1;
end

end
