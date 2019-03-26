function [imudata,imuparameter,imureference]=imudata_preprocessing(IMUdata,pos0,v0,Hz,eb,db,web,wdb)
%读取原始数据，并将数据设定为固定格式
% 输入：
%   IMUdata: imu采集数据,格式为：(Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z,Roll,Pitch,Yaw)
%   pos0:初始位置[lat;lon;hgt]输入单位度和米；	v0:初始速度m/s
%   Hz:采样频率;	
%     eb - gyro constant bias (deg/h)   
%     db - acc constant bias (ug)
%     web - angular random walk (deg/sqrt(h))
%     wdb - velocity random walk (ug/sqrt(Hz))

% 输出：
%   imudata:{gyr_x,gyr_y,gyr_z,acc_x,acc_y,acc_z,time} 单位为：弧度/s、m/s、s
%	imuparameter:imu原始数据参数，包括：
%       ts:采样周期；		Hz:采样频率;		L:数据长度(一般取整秒)
%       eb	db	web	wdb 最后全部转化为标准单位
%       pos0:初始位置[lat;lon;hgt]输入单位度和米,转化为弧度；	v0:初始速度m/s
%   imureference :若是IMUdata数据列数为9，则后面均为高精度的参考姿态数据，单位度



    imuparameter = [];
    imuparameter.Hz = Hz;
    imuparameter.ts = 1/Hz;
    L = fix(length(IMUdata)/Hz)*Hz;
    imuparameter.L = L;
    imuparameter.pos0 = pos0;
    imuparameter.v0 = v0;

    imuparameter.imuerr = imuerrset(eb, db, web, wdb); 
 
    imudata = zeros(L,6);
    imudata(1:L,1:3) = IMUdata(1:L,4:6);
    imudata(1:L,4:6) = IMUdata(1:L,1:3);
    
    [m,n] = size(IMUdata);
    if (n == 9)
        imureference = zeros(L,3);
        imureference = IMUdata(1:L,7:9);
    else
        imureference = [];
    end

