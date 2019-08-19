function data = dispData(f)
%DISPDATA   Useful information for DISPLAY at higher levels.
%   DATA = DISPDATA(F) extracts useful information from the given SINGFUN F and
%   the information DATA will be used by DISPLAY at higher levels. Currently, 
%   the only information it extracts is exponents.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% data{1}.name = 'exponents';
% data{1}.data = f.exponents;

data{1}.name = '  endpoint exponents';
exps = f.exponents;
data{1}.data = ['        ' '[' num2str(exps(1), '%2.2g') '      ' ...
            num2str(exps(2), '%2.2g') ']' '   '];

% More information can be appended to DATA:

end
