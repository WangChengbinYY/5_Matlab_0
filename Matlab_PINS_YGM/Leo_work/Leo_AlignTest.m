% 单子样、双子样，以及多子样的精度比对
%
%% 
%--------------------------------------------
% 初始化全局变量 glv，并获取imu采集数据 
% imu(MTI-G)数据格式为： (Acc_X，Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z,Roll,Pitch,Yaw)
%--------------------------------------------
    clear all;
    glvs
    [rpath, dpath, mytestflag] = psinsenvi();
%------------------读取mat格式的数据-----------------------------------  
% imudata:{gyr_x,gyr_y,gyr_z,acc_x,acc_y,acc_z,time} 单位为：弧度/s、m/s、s
    dpath = [dpath, '\IMUtest_short.mat'];
    load(dpath);
    Hz = 100;
    ts = 1/Hz;
    L = fix(length(IMUdata)/Hz)*Hz;
    imudata = zeros(L,3);
    imudata(1:L,1:3) = IMUdata(1:L,4:6);
    attrefer = zeros(L,4);
    attrefer(1:L,1:3) = IMUdata(1:L,7:9);    %参考姿态，单位 度
    attrefer(1:L,4) = [ts:ts:L*ts];
    
    
% 取1分钟数据平均，作为零偏，去除。
    imubias = mean(imudata(1:60/ts,1:3));
    imudata = imudata - imubias;    
    
    att0 =[0;0;0]*glv.deg;       
      
%%
%   第一种方法，直接单子样解算
    [qnb, att, Cnb] = attsyn(att0);  
    attfirst = zeros(L,4);
    attfirst(1,1:3) = att';
    attfirst(1,4) = ts;
    for i=2:L
        delta_phi = imudata(i-1,1:3)'*ts;
        Cbm_1bm = (eye(3)+askew(delta_phi));
        Cnb = Cnb*Cbm_1bm;
        Cnb = mnormlz(Cnb);                
        [qnb, att, Cnb] = attsyn(Cnb);
        attfirst(i,1:3) = att';
        attfirst(i,4) = attfirst(i-1,4) + ts;
    end
    attfirst(:,1:3) = attfirst(:,1:3)/glv.deg;
    attfirst(:,4) = attfirst(:,4)/60;
    
%%
%   第二种方法，双子样解算    
    [qnb, att, Cnb] = attsyn(att0);  
    nn = 2; nts = nn*ts;
    Len = fix(L/nn);
    attsecond = zeros(Len,4);
    ki = 1;
    for i = 1:nn:L-nn+1
        if i==1
            attsecond(1,1:3) = att';
            attsecond(1,4) = ts;
            ki = ki+1;
        else
            delta_phi1 = imudata(i-nn,1:3)'*ts;
            delta_phi2 = imudata(i-nn+1,1:3)'*ts;
            delta_phi = delta_phi1 + delta_phi2 +2*cros(delta_phi1,delta_phi2)/3;
            Cbm_1bm = rv2m(delta_phi);
            Cnb = Cnb*Cbm_1bm;            
            Cnb = mnormlz(Cnb);            
            [qnb, att, Cnb] = attsyn(Cnb);
            attsecond(ki,1:3) = att';  
            attsecond(ki,4) = attsecond(ki-1,4) + nts;
            ki = ki+1;
        end        
    end
    attsecond(:,1:3) = attsecond(:,1:3)/glv.deg;
    attsecond(:,4) = attsecond(:,4)/60;

