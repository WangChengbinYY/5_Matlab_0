function kf = kfUpdate_FtkfUpdate_FtPhiQ(kf,ins)
% 根据新的ins数据 更新kf的参数
% 2017.11.27 王成宾

 [wie,lat,RMh,RNh] = ... % just for short
     setvals(ins.eth.wie,ins.pos(1,1),ins.eth.RMh+ins.pos(3,1),ins.eth.RNh+ins.pos(3,1));
%先计算连续系统的Ft 
Ft = zeros(kf.n,kf.n);
Ft(1,2) = 2*wie*sin(lat);        Ft(2,1) = -Ft(1,2);
Ft(1,3) = -2*wie*cos(lat);       Ft(3,1) = -Ft(1,3);   
Ft(1,5) = -ins.eth.g;            Ft(2,4) = ins.eth.g;
Ft(4,2) = -1/RMh;                Ft(5,1) = 1/RNh;         Ft(6,1) = tan(lat)/RNh; %后面可以试验，看影响大不！！？？
Ft(4,5) = wie*sin(lat);          Ft(5,4) = -Ft(4,5);
Ft(4,6) = wie*cos(lat);          Ft(6,4) = -Ft(4,6);
Ft(1:3,10:12) = ins.Cnb;      Ft(4:6,7:9) = -ins.Cnb;      %是否取Cnb 或者 eye(3) 计算没有影响！！！！！！
kf.Ft = Ft;

%再计算 离散系统的
kf.Phikk_1 = kf.I + ins.nts*Ft + ins.nts^2*Ft^2/2;
M1 = kf.Qt;
M2 = Ft*M1 + (Ft*M1)';
kf.Qk = ins.nts*M1 + ins.nts^2*M2/2;

       