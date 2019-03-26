function kf = kfUpdate(kf,Zk)
% 根据新的观测数据Zk 进行Kalman滤波
% 2017.11.27 王成宾


     kf.Xkk_1 = kf.Phikk_1 * kf.Xk;
     kf.Pkk_1 = kf.Phikk_1 * kf.Pk * kf.Phikk_1' + kf.Qk;

     kf.Zk = Zk;
     kf.Kk = kf.Pkk_1*kf.Hk'*invbc(kf.Hk*kf.Pkk_1*kf.Hk'+kf.Rk);

     kf.XkZk = kf.Kk*(kf.Zk - kf.Hk*kf.Xkk_1);
     kf.Xk = kf.Xkk_1 + kf.XkZk;
     kf.Pk = (kf.I - kf.Kk*kf.Hk)*kf.Pkk_1*(kf.I - kf.Kk*kf.Hk)' + kf.Kk*kf.Rk*kf.Kk';




