function outPref = techPref(inPref)
%TECHPREF   Preference settings for CHEBTECH.
%   P = CHEBTECH.TECHPREF() returns a structure P with fields which contain the
%   default CHEBTECH preferences as field/value pairs.  This structure may be
%   passed to the CHEBTECH constructor.
%
%   P = CHEBTECH.TECHPREF(Q) does the same but replaces the default CHEBTECH
%   preferences with the values specified by the field/value pairs in the input
%   structure Q.
%
%   Note that no checking of either the input preference names or values takes
%   place and that preference names are case-sensitive.
%
%   ABSTRACT PREFERENCES REQUIRED OF ALL TECHS
%
%     eps          - Relative tolerance used in construction and subsequent
%      [2^-52]       operations.  See CHEBTECH.HAPPINESSCHECK for more details.
%
%     maxLength    - This preference is mapped to maxPoints (see below).
%
%     exactLength  - This preference is mapped to numPoints (see below).
%
%     extrapolate
%       true       - Function values at endpoints may be extrapolated from
%                    interior values rather than sampled.
%      [false]     - Do not extrapolate values at endpoints.
%
%     sampleTest
%      [true]      - Tests the function at one more arbitrary point to
%                    minimize the risk of missing signals between grid
%                    points.
%       false      - Do not test.
%
%   CHEBTECH-SPECIFIC PREFERENCES
%
%     gridType     - Type of Chebyshev grid on which the function is sampled.
%       1          - Use a grid of first-kind Chebyshev points.
%      [2]         - Use a grid of second-kind Chebyshev points.
%
%     minPoints    - Minimum number of points used by the constructor.  Should
%      [9]           be of the form 2^n + 1.  If not, it is rounded as such.
%
%     maxPoints    - Maximum number of points used by the constructor.
%      [2^16+1]
%
%     numPoints    - Fixed number of points used by constructor.  NaN allows
%      [NaN]         adaptive construction.
%
%     refinementFunction - Define function for refining sample values.
%      ['nested']        - Use the default process (nested evaluation).
%       'resampling'     - Every function value is computed afresh as the
%                          constructor tries grids of size 9, 17, 33, etc.
%       function_handle  - A user-defined refinement function.  See REFINE.m
%
%     happinessCheck     - Define function for testing happiness.
%      ['classic']       - Use the default process from Chebfun v4.
%       'strict'         - Strict tolerance for coefficients.
%       'loose'          - A looser tolerance for coefficients.
%       function_handle  - A user defined happiness. See HAPPINESSCHECK.m
%
% See also CHEBTECH, CHEBTECH1, CHEBTECH2

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

outPref.eps                = 2^-52;
outPref.gridType           = 2;
outPref.minPoints          = 9;
outPref.maxPoints          = 2^16 + 1;
outPref.numPoints          = NaN;
outPref.extrapolate        = false;
outPref.sampleTest         = true;
outPref.refinementFunction = 'nested';
outPref.happinessCheck     = 'classic';

if ( nargin == 1 )
    map.maxLength = 'maxPoints';
    map.exactLength = 'numPoints';
    outPref = chebfunpref.mergePrefs(outPref, inPref, map);
end

end
