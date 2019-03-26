function Qt = kfSetQt(kf,ins,imuerr)
% 依据ins参数，设定连续系统的 系统方差矩阵 Qt
% 2017.11.27 王成宾

%   w(t)=[加速度计常值零偏；陀螺常值零偏；后面都为0] 系统噪声
%   对应的方差强度矩阵，直接选取加速度计和陀螺对应的随机游走误差
wdb = imuerr.wdb;
web = imuerr.web;
Cnb = ins.Cnb;
Qt = zeros(kf.n,kf.n);
Qt(1,1) = Cnb(1,1)^2*wdb(1,1)^2 + Cnb(1,2)^2*wdb(2,1)^2 + Cnb(1,3)^2*wdb(3,1)^2;
Qt(1,2) = Cnb(1,1)*Cnb(2,1)*wdb(1,1)^2 + Cnb(1,2)*Cnb(2,2)*wdb(2,1)^2 + Cnb(1,3)*Cnb(2,3)*wdb(3,1)^2;
Qt(1,3) = Cnb(1,1)*Cnb(3,1)*wdb(1,1)^2 + Cnb(1,2)*Cnb(3,2)*wdb(2,1)^2 + Cnb(1,3)*Cnb(3,3)*wdb(3,1)^2;
Qt(2,1) = Qt(1,2);
Qt(2,2) = Cnb(2,1)^2*wdb(1,1)^2 + Cnb(2,2)^2*wdb(2,1)^2 + Cnb(2,3)^2*wdb(3,1)^2;
Qt(2,3) = Cnb(2,1)*Cnb(3,1)*wdb(1,1)^2 + Cnb(2,2)*Cnb(3,2)*wdb(2,1)^2 + Cnb(2,3)*Cnb(3,3)*wdb(3,1)^2;
Qt(3,1) = Qt(1,3);
Qt(3,2) = Qt(2,3);
Qt(3,3) = Cnb(3,1)^2*wdb(1,1)^2 + Cnb(3,2)^2*wdb(2,1)^2 + Cnb(3,3)^2*wdb(3,1)^2;

Qt(4,4) = Cnb(1,1)^2*web(1,1)^2 + Cnb(1,2)^2*web(2,1)^2 + Cnb(1,3)^2*web(3,1)^2;
Qt(4,5) = Cnb(1,1)*Cnb(2,1)*web(1,1)^2 + Cnb(1,2)*Cnb(2,2)*web(2,1)^2 + Cnb(1,3)*Cnb(2,3)*web(3,1)^2;
Qt(4,6) = Cnb(1,1)*Cnb(3,1)*web(1,1)^2 + Cnb(1,2)*Cnb(3,2)*web(2,1)^2 + Cnb(1,3)*Cnb(3,3)*web(3,1)^2;
Qt(5,4) = Qt(4,5);
Qt(5,5) = Cnb(2,1)^2*web(1,1)^2 + Cnb(2,2)^2*web(2,1)^2 + Cnb(2,3)^2*web(3,1)^2;
Qt(5,6) = Cnb(2,1)*Cnb(3,1)*web(1,1)^2 + Cnb(2,2)*Cnb(3,2)*web(2,1)^2 + Cnb(2,3)*Cnb(3,3)*web(3,1)^2;
Qt(6,4) = Qt(4,6);
Qt(6,5) = Qt(5,6);
Qt(6,6) = Cnb(3,1)^2*web(1,1)^2 + Cnb(3,2)^2*web(2,1)^2 + Cnb(3,3)^2*web(3,1)^2;