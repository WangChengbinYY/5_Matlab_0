function [Qt] = Leo_GetQtFrom(n,Cnb,web,wdb)
%计算新时刻的Qt
% 按照陀螺 加计的随机游走系数web wdb 来计算 方差强度阵q
% ！！！此处的web wdb 随机游走系数已经做过 单位换算了！！！
% 10/20/2017 Modified by Chengbin Wang, Tsinghua University
Qt = zeros(n);
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





