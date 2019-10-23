function kf = kfinit0(kf, nts)
% Always called by kfinit and initialize the remaining fields of kf.
    kf.nts = nts;
    [kf.m, kf.n] = size(kf.Hk);
    kf.I = eye(kf.n);
    if ~isfield(kf, 'xk'),  kf.xk = zeros(kf.n, 1);  end
    if ~isfield(kf, 'Qk'),  kf.Qk = kf.Qt*kf.nts;  end
    if ~isfield(kf, 'fading'),  kf.fading = 1;  end
    if ~isfield(kf, 'adaptive'),  kf.adaptive = 0;  end
    if kf.adaptive==1
        if ~isfield(kf, 'b'),  kf.b = 0.9;  end
        if ~isfield(kf, 'beta'),  kf.beta = 1;  end
        if ~isfield(kf, 'Rmin'),  kf.Rmin = 0.01*kf.Rk;  end
        if ~isfield(kf, 'Rmax'),  kf.Rmax = 100*kf.Rk;  end
        if ~isfield(kf, 'Qmin'),  kf.Qmin = 0.01*kf.Qk;  end
        if ~isfield(kf, 'Qmax'),  kf.Qmax = 100*kf.Qk;  end
    end
    if ~isfield(kf, 'xtau'),  kf.xtau = ones(size(kf.xk))*eps;   end
    if ~isfield(kf, 'T_fb'),  kf.T_fb = 1;   end
    kf.xfb = zeros(kf.n, 1);
%     kf.coef_fb = (1.0-exp(-kf.T_fb./kf.xtau));
%     kf.coef_fb = ar1coefs(kf.T_fb, kf.xtau);
    xtau = kf.xtau;
    xtau(kf.xtau<kf.T_fb) = kf.T_fb;  kf.coef_fb = kf.T_fb./xtau;  %2015-2-22