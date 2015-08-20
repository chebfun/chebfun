function v_cheb = legcoeffs2chebvals(c_leg, varargin)

c_cheb = leg2cheb(c_leg, varargin{:});
v_cheb = chebvals2chebcoeffs(c_cheb);

end