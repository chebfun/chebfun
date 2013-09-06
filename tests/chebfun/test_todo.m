function pass = test_todo(pref)

% Ensure we've written all the test files for this branch.

files = {'test_airy' % Done
'test_angle' % Requires ATAN2()
% 'test_area' % Not rquired
'test_besselh'
'test_besseli'
'test_besselj' % Requires POWER()
'test_besselk'
'test_bessely'
'test_cov' % Requires MEAN() which requires SUM
'test_cumprod' % Not required?
'test_ellipj'
'test_ellipke'
'test_end' % Requires SUBSREF()
% 'test_erfcinv' % Done test_erfX
% 'test_erfc' % Done test_erfX
% 'test_erfcx' % Done test_erfX
% 'test_erfinv' % Done test_erfX
% 'test_erf' % Done test_erfX
% 'test_fill' % Not rquired
'test_hypot' % Requires addBreaksAtRoots()
'test_interp1'
% 'test_legpoly' % Not required
'test_mean' % requires SUM()
'test_pchip'
% 'test_poly' % Not required
'test_prod' % Not required?
'test_residue'
'test_spline'
'test_std' % Requires VAR()
'test_var'}; % Requires MEAN() which requires SUM()

for k = 1:numel(files)
    pass = exist(files{k}, 'file');
    if ( ~pass )
        break
    end
end

% [TODO]: Delete this file.