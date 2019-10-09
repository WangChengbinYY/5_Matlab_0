% SINS/GPS intergrated navigation simulation unsing kalman filter.
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS, test_SINS_GPS_186, test_SINS_GPS_193.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/06/2011
clear;

glvs
psinstypedef(153);
trj = trjfile('trj10ms.mat');
% initial settings
[nn, ts, nts] = nnts(1, trj.ts);

% 读取原始IMU数据，这里添加随机噪声一次后，保存数据
% load('trj10ms_imu.mat');
% load('trj10ms_imuerr.mat');
% imuerr = imuerrset(0.03, 100, 0.001, 5);
imuerr = imuerrset(2.5, 13, 0.15, 3700);   %ADIS16467
imu = imuadderr(trj.imu, imuerr);

% 设定起始位置的初始误差 姿态单位 前面两个是角秒，航向是 角分；速度是m/s，位置是 m
% load('trj10ms_davp0.mat');
% load('trj10ms_avp0_err.mat');
davp0 = avpseterr([0.5*3600;-0.5*3600;60], 0.5, [3;3;5]);
avp0_err = avpadderr(trj.avp0,davp0);
avp0 = trj.avp0;
% 不含误差的精确数值为
% load('trj10ms_avp_theory.mat');
avp_theory = trj.avp;

% 利用精确地位置信息，构建有误差的GPS数据
% load('trj10ms_GPS.mat');
GPS = zeros(length(avp_theory),4);
GPS(:,4) = avp_theory(:,10);
for i = 1:length(avp_theory)
    GPS(i,1:3) = avp_theory(i,7:9)' + davp0(7:9).*randn(3,1); 
end

ins = insinit(avp0_err, ts);  ins.nts = nts;

% KF filter
kf = kfinit(ins, davp0, imuerr);
len = length(imu); 

avp_result_LC = zeros(len,10);

for i=1:len
    wvm = imu(i,1:6);  t = imu(i,end);
    ins = insupdate(ins, wvm);
    avp_result_LC(i,:) = [ins.avp', t];
    kf.Phikk_1 = kffk(ins);
    kf = kfupdate(kf);
    if mod(i,100)==0
        kf = kfupdate(kf, ins.pos-GPS(i,1:3)', 'M');
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');
        avp_result_LC(i,:) = [ins.avp', t];
    end
end

% Result_Zk = zeros(fix(len/100),6);
% j=1;
% 
% kf.Phikk_1 = eye(15);
% for i=1:len
%     wvm = imu(i,1:6);  t = imu(i,end);
%     ins = insupdate(ins, wvm);
%     avp_result_LC(i,:) = [ins.avp', t];
%     Phikk_1 = kffk(ins);
%     kf.Phikk_1 = Phikk_1*kf.Phikk_1;
%     if mod(i,100)==0
%         Zk=ins.pos-GPS(i,1:3)';
%         kf = kfupdate(kf, Zk, 'B');
%         [kf, ins] = kffeedback(kf, ins, 1, 'vp');
%         avp_result_LC(i,:) = [ins.avp', t];
%         kf.Phikk_1 = eye(15);
%         
%         Result_Zk(j,1:3) = Zk';
%         Result_Zk(j,4:6) = GPS(i,1:3);
%         j=j+1;
%     end
% end


% show results
avperr = avpcmp(avp_result_LC, avp_theory);
insplot(avp_result_LC);
inserrplot(avperr);
% kfplot(xkpk, avperr, imuerr);

