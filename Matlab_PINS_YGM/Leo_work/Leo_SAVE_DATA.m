% 生成用于初始对准Kalman滤波的数据
% 10/20/2017 Modified by Chengbin Wang, Tsinghua University
clear all;
glvs
ts = 0.1;   % sampling interval
T = 1000;
avp0 = avpset([3;2;30], [0;0;0], [30;108;380]);
imuerr = imuerrset(1, 100, 0.1, 10);
imu = imustatic(avp0, ts, T, imuerr);   % IMU simulation

[rpath, dpath, mytestflag] = psinsenvi();
dpath = [dpath, '\TestData_Align.mat'];
parameter = struct('ts',ts,'T',T,'attitude0',avp0(1:3)/glv.deg,'velocity0',...
    avp0(4:6),'position0',[avp0(7:8)/glv.deg;avp0(9)],'eb',imuerr.eb/glv.dph,...
    'db',imuerr.db/glv.ug,'web',imuerr.web/glv.dpsh,'wdb',imuerr.wdb/glv.ugpsHz);

describe = [sprintf('IMU仿真数据参数为：\n'),...
    sprintf('   仿真时间：采样 %.3fs，周期 %.3fs；\n',parameter.ts,parameter.T),...
    sprintf('   真实姿态：俯仰 %.3f度，横滚 %.3f度，航向 %.3f度；\n',avp0(1:3)/glv.deg),...
    sprintf('   初始速度：东向 %.3fm/s，北向 %.3fm/s，天向 %.3fm/s；\n',avp0(4:6)),...
    sprintf('   初始位置：纬度 %.3f度，经度 %.3f度，高程 %.3fm；\n',avp0(7:8)/glv.deg,avp0(9)),...
    sprintf('   陀螺误差：零偏 %.3f度/h，角度随机游走 %.3fdeg/sqrt(h)；\n',imuerr.eb(1)/glv.dph,imuerr.web(1)/glv.dpsh),...
    sprintf('   加计误差：零偏 %.3fug，速率随机游走 %.3fug/sqrt(Hz)。\n',imuerr.db(1)/glv.ug,imuerr.wdb(1)/glv.ugpsHz)];

save(dpath,'imu','avp0','imuerr','describe','ts','T');
%imu - gyro & acc incremental outputs
%imu(1:3) gyro 单位弧度 imu(4:6) acc 单位m/s

% clear all;

% msgbox('The simulation of IMU is done!');      


 [attsb, qnb] = alignsb(imu, avp0(7:9));
  attsb/glv.deg
  
davp = avpseterr([-30;30;30], [0.01;0.01;0.01]*0, [1;1;1]*0);
avp = avpadderr(avp0, davp);
  
 wvn = [0.01;0.01;0.01];
 phi0 = [0.5;0.5;1]*glv.deg;

[att0v, attkv, xkpkv] = alignvn(imu, [3.4;2.6;30.8]*glv.deg, avp0(7:9), phi0, imuerr, wvn);
att0v/glv.deg
[att0v, attkv, xkpkv] = alignvn1(imu,[3.4;2.6;30.8]*glv.deg, avp0(7:9), phi0, imuerr, wvn);
att0v/glv.deg

