% ѧϰ�Ϲ��� ���ٶȷ��� kalman�˲�


%--------------------------------------------
% ��ʼ��ȫ�ֱ��� glv������ȡ�������� 
% imu ���� ���� �Ӽ� ʱ�� ����
% parameter ���� ���ַ������
%--------------------------------------------
    clear all;
    glvs
    [rpath, dpath, mytestflag] = psinsenvi();
    dpath = [dpath, '\TestData_Align.mat'];
    load(dpath);
    %disp(describe);
    
[attsb, qnb] = alignsb(imu, avp0(7:9));
wvn = [0.01;0.01;0.01];
phi0 = [0.5;0.5;5]*glv.deg;
[att0v, attkv, xkpkv] = alignvn1(imu, attsb, avp0(7:9), phi0, imuerr, wvn);
attsb/glv.deg
att0v/glv.deg
attkv(5000,:)/glv.deg


    
    
    
    
    
    
    
    
    
    
    