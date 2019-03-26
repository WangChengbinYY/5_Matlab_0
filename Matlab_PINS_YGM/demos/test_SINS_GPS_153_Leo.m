% SINS/GPS intergrated navigation simulation unsing kalman filter.
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS, test_SINS_GPS_186, test_SINS_GPS_193.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/06/2011
clear variables;

glvs
psinstypedef(153);
% trj = trjfile('trj10ms.mat');
% initial settings
% [nn, ts, nts] = nnts(1, trj.ts);
nn = 1;
ts = 1/200;
nts = nn*ts;
imuerr = imuerrset(0.01, 10.2246072, 0.0022, 1.278075906068187,0.005,4,25,4);   %改成采集的数据参数
% imuerr = imuerrset(0.01, 10, 5, 1.25);
% imu = imuadderr(trj.imu, imuerr);
% davp0 = avpseterr([30;-30;20], 0.1, [1;1;3]);
%初始位置，含有误差 [姿态；速度；位置]'
% 初始姿态误差，取 水平 0.5度 航向 1度 试试看吧
% 速度误差 取 0.1m/s2
davp0 = [0.5*pi/180;0.5*pi/180;-3.1259987;0.1;0.1;0.1;0.530699173342059;1.994351183398699;15.359326045134754];
davp0Error = [0.5*pi/180;0.5*pi/180;1*pi/180;0.1;0.1;0.1;5/6378137;5/6378137;5];


ins = insinit(davp0, ts);  ins.nts = nts;
ins.tauG = [4*3600;4*3600;4*3600];
ins.tauA = [4*3600;4*3600;4*3600];
ins.eb = imuerr.eb;
ins.db = imuerr.db;
% KF filter
kf = kfinit(ins, davp0Error, imuerr);

%读取IMU 和 GPS数据
load('E:\n_WorkSpace\Matlab\0_Data\INSGPS组合使用数据\INS_GPS');
imu(:,1:6) = imudata(:,2:7)*ts;
% imu(:,3) = (-1)*imu(:,3);
imu(:,7) = imudata(:,1);
NUM_GNSS = 1;


len = length(imu); 
xkpk = prealloc(fix(len*ts), 2*kf.n+1);
avp = prealloc(len, 10);
bias = prealloc(fix(len*ts),7);
%测试用
test_pos = zeros(len,10);
test_xk  = zeros(fix(len*ts),10);

% timebar(nn, len, '15-state SINS/GPS Simulation.'); 
ki = 1;
% tbstep = floor(len/nn/100); tbi = timebar(1, 99, 'SINS/GPS Simulation.');
profile on
for k=1:nn:len-nn+1
    k1 = k+nn-1;  
    wvm = imu(k:k1,1:6);  t = imu(k1,end);
    ins = insupdate(ins, wvm);
    kf.Phikk_1 = kffk(ins);
    kf = kfupdate(kf);  
    
    tmp_time = abs(imu(k,7)-GNSS(NUM_GNSS,1));
    if( tmp_time < 0.004)
        posGPS(1:2,1) = GNSS(NUM_GNSS,2:3)'*pi/180;
        posGPS(3,1) = GNSS(NUM_GNSS,4);
        NUM_GNSS = NUM_GNSS +1;
        DeltaPos =  ins.pos-posGPS;
        kf = kfupdate(kf,DeltaPos, 'B');
        
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t]'; 
        [kf, ins] = kffeedback(kf, ins, 1, 'avped');       
        bias(ki,1:3) = ins.eb';
        bias(ki,4:6) = ins.db';
        bias(ki,7) = t;
        
        ki = ki+1;
    end
    
    avp(k1,:) = [ins.avp',t];     
    

end
profile viewer
% show results
% avperr = avpcmp(avp, trj.avp);
insplot(avp);
% inserrplot(avperr);
% kfplot(xkpk, imuerr);

