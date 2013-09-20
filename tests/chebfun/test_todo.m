function pass = test_todo(pref)

% Ensure we've written all the test files for this branch.

files = {
%     'test_airy' % DONE
% 'test_angle' % Not required.
% 'test_area' % Not rquired
% 'test_besselh'
% 'test_besseli' % Not required.
'test_besselj' % Requires POWER()
% 'test_besselk'
% 'test_bessely'
% 'test_cumprod' % Not required?
% 'test_ellipj'
% 'test_ellipke'
% 'test_end' % DONE.
% 'test_erfcinv' % DONE test_erfX
% 'test_erfc' % DONE test_erfX
% 'test_erfcx' % DONE test_erfX
% 'test_erfinv' % DONE test_erfX
% 'test_erf' % DONE test_erfX
% 'test_fill' % Not rquired
% 'test_hypot' % DONE.
% 'test_interp1' % DONE.
% 'test_legpoly' % Not required
% 'test_mean' % DONE.
% 'test_pchip' % DONE.
% 'test_poly' % Not required
% 'test_prod' % Not required?
% 'test_residue' % DONE
% 'test_spline' % DONE
% 'test_cov'
% 'test_std' % Requires VAR()
% 'test_var' % Requires MEAN() which requires SUM()
};

for k = 1:numel(files)
    pass = exist(files{k}, 'file');
    if ( ~pass )
        break
    end
end

% [TODO]: Delete this file.
