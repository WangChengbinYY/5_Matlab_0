function [imudata,imuparameter,imureference]=imudata_preprocessing(IMUdata,pos0,v0,Hz,eb,db,web,wdb)
%��ȡԭʼ���ݣ����������趨Ϊ�̶���ʽ
% ���룺
%   IMUdata: imu�ɼ�����,��ʽΪ��(Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z,Roll,Pitch,Yaw)
%   pos0:��ʼλ��[lat;lon;hgt]���뵥λ�Ⱥ��ף�	v0:��ʼ�ٶ�m/s
%   Hz:����Ƶ��;	
%     eb - gyro constant bias (deg/h)   
%     db - acc constant bias (ug)
%     web - angular random walk (deg/sqrt(h))
%     wdb - velocity random walk (ug/sqrt(Hz))

% �����
%   imudata:{gyr_x,gyr_y,gyr_z,acc_x,acc_y,acc_z,time} ��λΪ������/s��m/s��s
%	imuparameter:imuԭʼ���ݲ�����������
%       ts:�������ڣ�		Hz:����Ƶ��;		L:���ݳ���(һ��ȡ����)
%       eb	db	web	wdb ���ȫ��ת��Ϊ��׼��λ
%       pos0:��ʼλ��[lat;lon;hgt]���뵥λ�Ⱥ���,ת��Ϊ���ȣ�	v0:��ʼ�ٶ�m/s
%   imureference :����IMUdata��������Ϊ9��������Ϊ�߾��ȵĲο���̬���ݣ���λ��



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

