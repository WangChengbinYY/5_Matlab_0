% SINS pure inertial navigation simulation. Please run 
% 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS_GPS_153, test_SINS_static.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/06/2011

% 200Hz 四子样算法
clear all;
profile on
glvs

load('E:\2_WorkSpace_Leo\Matlab\Matlab_SINS\data\imudata_pos80');
Hz = 200;
%数据是 角速度和加速度————改成增量形式
imu(:,1:6) = imudata(586800:1153371,2:7)/Hz;
imu(:,7) = imudata(586800:1153371,1);

% trj = trjfile('trj10ms.mat');
[nn, ts, nts] = nnts(2, 1/Hz);
% imuerr = imuerrset(0.01, 100, 0.001, 10);
% imu = imuadderr(trj.imu, imuerr);
avp00 =zeros(9,1);
avp00(3)=-3.126;
avp00(7)=30.4068572011000*pi/180.0;
avp00(8)=114.267934323500*pi/180.0;
avp00(9)=20.3593264287000;
avp = inspure(imu, avp00);  % pure inertial navigation

