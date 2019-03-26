function kf = kfinit_wcb(n,m)
% 静基座对准，速度反馈Kalman参数 初始化
% 2017.11.27 王成宾

kf = [];
kf.n = n; % 状态方程维数
kf.m = m; % 观测方程维数
kf.I = eye(n);  %单位矩阵
kf.Xk_1 = zeros(n,1);
kf.Pkk_1 = zeros(n,n);
kf.Xk = zeros(n,1);
kf.Pk = zeros(n,n);
kf.Zk = zeros(m,1);
kf.Hk = zeros(m,n);
kf.Kk = zeros(n,m);
kf.XkZk = zeros(n,1);  %记录中间的 测量Zk对Xk的贡献 Kk(Zk-HkXkk_1);

kf.Ft = zeros(n,n);     % 连续系统的状态矩阵
kf.Qt = zeros(kf.n);    % 连续系统的状态方程系统噪声方差阵
kf.Phikk_1 = zeros(n,n);   % 离散化的状态矩阵
kf.Qk = zeros(n,n);
kf.Rk = zeros(m,m);



