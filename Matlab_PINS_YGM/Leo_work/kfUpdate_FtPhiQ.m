function kf = kfUpdate_FtkfUpdate_FtPhiQ(kf,ins)
% �����µ�ins���� ����kf�Ĳ���
% 2017.11.27 ���ɱ�

 [wie,lat,RMh,RNh] = ... % just for short
     setvals(ins.eth.wie,ins.pos(1,1),ins.eth.RMh+ins.pos(3,1),ins.eth.RNh+ins.pos(3,1));
%�ȼ�������ϵͳ��Ft 
Ft = zeros(kf.n,kf.n);
Ft(1,2) = 2*wie*sin(lat);        Ft(2,1) = -Ft(1,2);
Ft(1,3) = -2*wie*cos(lat);       Ft(3,1) = -Ft(1,3);   
Ft(1,5) = -ins.eth.g;            Ft(2,4) = ins.eth.g;
Ft(4,2) = -1/RMh;                Ft(5,1) = 1/RNh;         Ft(6,1) = tan(lat)/RNh; %����������飬��Ӱ��󲻣�������
Ft(4,5) = wie*sin(lat);          Ft(5,4) = -Ft(4,5);
Ft(4,6) = wie*cos(lat);          Ft(6,4) = -Ft(4,6);
Ft(1:3,10:12) = ins.Cnb;      Ft(4:6,7:9) = -ins.Cnb;      %�Ƿ�ȡCnb ���� eye(3) ����û��Ӱ�죡����������
kf.Ft = Ft;

%�ټ��� ��ɢϵͳ��
kf.Phikk_1 = kf.I + ins.nts*Ft + ins.nts^2*Ft^2/2;
M1 = kf.Qt;
M2 = Ft*M1 + (Ft*M1)';
kf.Qk = ins.nts*M1 + ins.nts^2*M2/2;

       