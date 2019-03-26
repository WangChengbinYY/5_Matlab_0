% 11/22/2017 Modified by Chengbin Wang, Tsinghua University
%% 
%--------------------------------------------
% 初始化全局变量 glv，并获取imu采集数据 
% imu(MTI-G)数据格式为： (Acc_X，Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z,Roll,Pitch,Yaw)
%--------------------------------------------
    clear all;
    glvs
    [rpath, dpath, mytestflag] = psinsenvi();
%------------------读取mat格式的数据-----------------------------------    
     dpath = [dpath, '\IMUdataNew'];
     load(dpath);
%%  初始化        
    % 设置IMU器件的误差参数，并进行单位转化   
    eb = 10; %     eb - gyro constant bias (deg/h)     
    db = 40; %     db - acc constant bias (ug)
    web = 0.6; %     web - angular random walk (deg/sqrt(h))
    wdb = 100;%     wdb - velocity random walk (ug/sqrt(Hz))
 
    % 设置初始位置  position [lat;lon;hgt] 和 初始 速度
    pos0 = [39.90*glv.deg;116.39*glv.deg;0];
    v0 = [0;0;0];
    Hz = 100; %采样频率
    % 初始化原始数据
    %   imudata:{gyr_x,gyr_y,gyr_z,acc_x,acc_y,acc_z,time} 单位为：弧度/s、m/s、s
    [imudata,imuparameter,imureference]=imudata_preprocessing(IMUdataNew,pos0,v0,Hz,eb,db,web,wdb);
    % 清理中间变量
%     clear eb db web wdb Hz v0 pos0 IMUdata;
    
%%       
%------------------------------------------
% 第一步 ：先利用1分钟的陀螺和加计数据进行 静基座 粗对准 获取初始姿态 attsb
%------------------------------------------
     [attsb, qnb] = alignsb(imudata,imuparameter.pos0);
     disp('静基座1分钟粗对准结果：');   
     attsb/glv.deg

     %C
%      ts = 0.01;
%     phi = [1;1;5]*glv.deg;   %假设的大失准角预估
%     wvn = [0.01;0.01;0.01];
%     [att0v, attkv, xkpkv] = alignvn(imudata, qnb, pos0, phi, imuparameter.imuerr, wvn,ts);
%     disp('速度反馈精对准结果：');   
%     att0v/glv.deg  
     
     
%%       
%------------------------------------------
% 第二步 ：设置速度反馈Kalman滤波参数的初始化
%------------------------------------------
    attsb = IMUdataNew(1,7:9)'*glv.deg;
    nn = 2; %设置姿态更新的子样数，这里先采用简单的 二子样算法,nts=nn*ts,表示INS的解算周期
    %设置惯性导航解算数据结构体，并初始化
    ins = insinit_wcb(attsb,imuparameter.v0,imuparameter.pos0,imuparameter.ts,nn);
    
    %设置速度反馈Kalman滤波参数的结构体 状态方程12维,观测量3维
    % Xt=[速度误差，姿态误差，陀螺误差，加速度计误差]
    % Zt=[速度]
    n = 12; m = 3;
    kf = kfinit_wcb(n,m);
    
    % Set X0 P0  Hk0 Rk0;  Ft0 Qt0
    kf.Xk = [[0.01;0.01;0.01];[0.5;0.5;5]*glv.deg;[0;0;0];[0;0;0]]; % x0
                % 初始速度误差    姿态误差       陀螺零偏     加计零偏
    kf.Pk = diag([[0.01;0.01;0.01];[0.5;0.5;5]*glv.deg;imuparameter.imuerr.eb;imuparameter.imuerr.db])^2;    % P0
    kf.Hk(1:3,1:3) = eye(3); 
    kf.Rk = diag([0.01;0.01;0.01])^2;   %观测噪声方差 因为是静基座，设定为0.01m/s2    
    %依据ins参数计算kf的Qt、Ft,再计算Phikk_1、Qk
    %kf.Qt = kfSetQt(kf,ins,imuparameter.imuerr);  %此处简化了，后面不需要更新了，原理，见书上
    kf.Qt = diag([imuparameter.imuerr.wdb;imuparameter.imuerr.web;[0;0;0];[0;0;0]])^2;
                %加速度计随机游走       陀螺随机游走
    kf = kfUpdate_FtPhiQ(kf,ins);          
    
%%    
%------------------------------------------
% 第三步 ：开始循环计算惯导解算以及Kalman滤波
%------------------------------------------    
% 设定所要存储的中间数据
    len = fix(imuparameter.L/ins.nn)*ins.nn;
    avp = zeros(len,10);     %存储中间计算的 姿态、速度、位置、时间
    ki = 1;
%=============================================
% 记录中间变量
    temp_vn_ins = zeros(len,3);     %C 记录惯导解算的速度
    temp_vn_kf = zeros(len,3);      %C 计算kf滤波的速度误差
    temp_att_ins = zeros(len,3);     %C 记录惯导解算的姿态
    temp_att_kf = zeros(len,3);      %C 计算kf滤波的姿态误差
    temp_Xk = zeros(len,12);        %C 
    temp_Pk = zeros(len,12);
    temp_Xkk_1 = zeros(len,12);
    temp_XkZk = zeros(len,12);      %C 测量Zk对Xk的贡献 Kk(Zk-HkXkk_1);
%=============================================
    for i=1:nn:len-nn+1
       %先更新惯导解算
       %先用两子样算法更新姿态
        if i==1
            avp(ki,1:3) = ins.att'/glv.deg;
            avp(ki,4:6) = ins.vn';    avp(ki,7:9) = ins.pos';
            avp(ki,10) = imuparameter.ts;
            temp_vn_ins(ki,:) = ins.vn';    %C
            ki = ki+1;            
        else
            %依据上一时刻的姿态和加速度计 计算此刻的速度、位置
            avp(ki,10) = avp(ki-1,10) + ins.nts;  %当前时刻的时间
            ins.eth = earth(ins.pos, ins.vn);   %更新出 wnie,wnen,gn
            %-----------kf部分参数更新(利用tm_1时刻的参数更新Phikk_1、Qk)-- 
            %kf.Qt = kfSetQt(kf,ins,imuparameter.imuerr); %已简化，不需要更新
            kf = kfUpdate_FtPhiQ(kf,ins);
            
            %-----------速度更新--------------------------------------
            delta_fb = (imudata(i-nn,4:6)+imudata(i-nn+1,4:6))'/2;
            delta_vn = ins.Cnb * delta_fb + cros(2*ins.eth.wnie+ins.eth.wnen,ins.vn)+...
                ins.eth.gn;
            ins.vn = ins.vn + delta_vn*ins.nts; 
            temp_vn_ins(ki,:) = ins.vn';    %C
            
            %C 
            if norm(delta_vn)>10
                temp_bug = 0;
            end
            
            %-----------姿态更新--------------------------------------
            delta_phi1 = imudata(i-nn,1:3)'*ins.ts;
            delta_phi2 = imudata(i-nn+1,1:3)'*ins.ts;
            delta_phi = delta_phi1 + delta_phi2 +2*cros(delta_phi1,delta_phi2)/3;
            Cbm_1bm = rv2m(delta_phi);
            ins.Cnb = ins.Cnb*Cbm_1bm;            
            ins.Cnb = mnormlz(ins.Cnb);        %Cnb单位正交化    
            [ins.qnb, ins.att, ins.Cnb] = attsyn(ins.Cnb);
            temp_att_ins(ki,:) = ins.att';     %C
            
            %-----------以更新后的速度为误差观测量，进行Kalman滤波
            kf = kfUpdate(kf,ins.vn);
            temp_vn_kf(ki,:) = kf.Xk(1:3,1)';   %C
            temp_att_kf(ki,:) = kf.Xk(4:6,1)';  %C
            %C 并记录中间的过程变量
            temp_Xk(ki,:) = kf.Xk';        %C 
            temp_Pk(ki,:) = diag(kf.Pk)';        %C 
            temp_Xkk_1(ki,:) = kf.Xkk_1';        %C 
            temp_XkZk(ki,:) = kf.XkZk';      %C 测量Zk对Xk的贡献 Kk(Zk-HkXkk_1);
            %----------利用Kalman滤波结果，更新姿态及速度信息
            ins.qnb = qdelphi(ins.qnb, 0.1*kf.Xk(4:6,1));
            [ins.qnb, ins.att, ins.Cnb] = attsyn(ins.qnb);
            kf.Xk(4:6,1) = 0.9*kf.Xk(4:6,1);
            ins.vn = ins.vn - 0.1*kf.Xk(1:3,1);
            kf.Xk(1:3,1) = 0.9*kf.Xk(1:3,1); 
            avp(ki,1:3) = ins.att'/glv.deg;
            avp(ki,4:6) = ins.vn'; 
            
            ki = ki+1;
            
            %C  调试设定的断点
            if ki == 10240
                temp_a = 0;
            end
            
        end       
    end
    
    temp_att_ins = temp_att_ins/glv.deg;



%     avp0 = zeros(1,10);
%     avp0(1:3)=attsb';
%     avp0(7:9)=pos0';
%     avp0(10)=ts;
% 
%     ins = insinit(avp0,ts);
%     nn=2; len=L;
%     avp = zeros(len,10); kk = 1;
%     for k= 1:nn:len-nn+1
%         k1 = k+nn-1;
%         wvm = IMUdata(k:k1,1:6);   t = k*ts;
%         ins = insupdate(ins, wvm);
%         avp(kk,:) = [q2att(ins.qnb);ins.vn;ins.pos;t]';
%         kk = kk+1;        
%     end
%     
%     
%     
% 
% 
% 
% %% 建立Kalman滤波方程变量的 初始化
%     nn = 2;  %  设定 Kalman滤波周期内的采样数 也就是 n子样计算法
%     
%     
%     
%     
% %% vn-meas. Kalman filter
% %-----------------------------------------
% % 利用严恭敏的 速度反馈kalman滤波 计算结果为
% %-----------------------------------------
% phi = [1;1;5]*glv.deg;   %假设的大失准角预估
% wvn = [0.01;0.01;0.01];
% [att0v, attkv, xkpkv] = alignvn(IMUdata, qnb, pos0, phi, imuerr, wvn,ts);
%     disp('速度反馈精对准结果：');   
%     att0v/glv.deg    
% 
%     
% %% gyro-compass method
% ctl0 = [20; 30]; ctl1 = [50; 300];
% [att0c, attkc] = aligncmps(IMUdata, qnb, pos0, ctl0, ctl1,ts);    
%     disp('速度反馈精对准结果：');   
%     att0c/glv.deg   
%     
% %%
% %------------------------------------------
% % 在静基座粗对准基础上开始 速度反馈的Kalman滤波精对准（不包含圆锥 划船补偿）
% %------------------------------------------
% 
% %% ----惯导解算初始化-----
% 
%     nn = 4; nts = nn*ts;                        %没隔nn个采样周期进行一次Kalman滤波
%     len = fix(length(imu)/nn)*nn;          %数据长度为len
%     [IMU_att,IMU_ven] = prealloc(len, 3, 3);    %惯导解算出来的姿态和速度数据 姿态单位为弧度
%     phi0 = [0.5;0.5;2]*glv.deg;                     %设定初始失准角
%     IMU_att(1,:) = attsb'; 
%     wvn = [0.01;0.01;0.01];                  %粗 设定初始速度    
%     IMU_ven(1,:) = wvn';
%     Cnb = a2mat(attsb); qnb = a2qua(attsb);
% 
% %% ----Kalman滤波初始化-----
%     eth = earth(avp0(7:9));                 %获取地球相关计算参数
%     v_error = zeros(3,1);  p_error = zeros(3,1); %速度、位置误差
% %第一步 进行陀螺、加计参数设定
%     % imuerr;
% %第二步 进行Kalmn滤波参数初始设定 包括Ft0,Qt0, Rt0, Pt0, Ht0及初始姿态
%     kf = [];
%     kf.nts = nts;
%     kf.Cnb = a2mat(attsb(1:3));     %利用粗对准结果
%     kf.qnb = a2qua(attsb(1:3));
%     kf.n = 12; kf.m = 3;            %n 状态量个数，m 观测量个数
%     kf.I = eye(kf.n);
%     kf.Xk = zeros(kf.n,1); kf.Zk = zeros(kf.m,1);
%     kf.Ft = zeros(kf.n); kf.Hk = zeros(kf.m,kf.n);
%     kf.Qk = zeros(kf.n); kf.Qt = zeros(kf.n); 
%     kf.Rk = zeros(kf.m);  kf.Pk = zeros(kf.n); 
%     kf.Pkk_1 = zeros(kf.n); kf.Xkk_1 = zeros(kf.n,1); 
%     kf.Kk = zeros(kf.n);  kf.Phikk_1 = zeros(kf.n);
% 
% %     kf.Xk = [wvn; phi0; [0;0;0]; [0;0;0]];   %系统初始状态 全设为0
%     kf.Pk = diag([[1;1;1]; phi0; imuerr.eb; imuerr.db])^2;      %初始 P0
%     kf.Rk = diag(IMU_ven(1,1))^2;
% 
%     [xkpk,xkpksb,xkpksb1] = prealloc(fix(len/nn), 2*kf.n, 2*kf.n, 2*kf.n);  %将中间系统状态 数据记录下来
%     [IMU_kfattk,IMU_kfvnk] = prealloc(fix(len/nn), 3, 3); 
% 
%     %先计算一下 Ft中的不变量
%     [wie,lat,RMh,RNh] = ... % just for short
%     setvals(glv.wie,avp0(7,1),eth.RMh+avp0(9,1),eth.RNh+avp0(9,1)); 
%     kf.Ft(1,2) = 2*wie*sin(lat);        kf.Ft(2,1) = -kf.Ft(1,2);
%     kf.Ft(1,3) = -2*wie*cos(lat);       kf.Ft(3,1) = -kf.Ft(1,3);   
%     kf.Ft(1,5) = -eth.g;                kf.Ft(2,4) = eth.g;
%     kf.Ft(4,2) = -1/RMh;                kf.Ft(5,1) = 1/RNh;         kf.Ft(6,1) = tan(lat)/RNh; %后面可以试验，看影响大不！！？？
%     kf.Ft(4,5) = wie*sin(lat);          kf.Ft(5,4) = -kf.Ft(4,5);
%     kf.Ft(4,6) = wie*cos(lat);          kf.Ft(6,4) = -kf.Ft(4,6);
%     kf.Ft(1:3,10:12) = kf.Cnb;          kf.Ft(4:6,7:9) = -kf.Cnb;      %是否取Cnb 或者 eye(3) 计算没有影响！！！！！！
% %     kf.Ft(1:3,10:12) = eye(3);          kf.Ft(4:6,7:9) = -eye(3);
%     kf.Hk(:,1:3) = eye(3);
%     
%      [kfZk_HkXkk_1,IMU_Xk_add] = prealloc(fix(len/nn), 6, 6); 
% %% ----开始循环惯导解算，每隔nn个采样周期 加入一次Kalman滤波-----
%     dven = zeros(3,1);  delphi = zeros(3,1);  mm = 1;
%     for i=1:1:len
%         %% 判断若是满足 nn 间隔，则加入 Kalman滤波  
%         if mod(i,nn) == 0
%              %计算新时刻的Ft
%               kf.Ft(1:3,10:12) = kf.Cnb; kf.Ft(4:6,7:9) = -kf.Cnb;
%              %计算新时刻的Qt
%              % 按照陀螺 加计的随机游走系数web wdb 来计算 方差强度阵q
%              % ！！！此处的web wdb 随机游走系数已经做过 单位换算了！！！,并且直接和系统状态方程相关，不通用
%              kf.Qt = Leo_GetQtFrom(kf.n,Cnb,imuerr.web,imuerr.wdb);  %后面做实验 直接用固定的，是否有变化             
%              [kf.Phikk_1, kf.Qk] = kfc2d(kf.Ft, kf.Qt, ts*nn, 3);  
%             
%              kf.Xkk_1 = kf.Phikk_1 * kf.Xk;
%              kf.Pkk_1 = kf.Phikk_1 * kf.Pk * kf.Phikk_1' + kf.Qk;
% %                 kf.Pkk_1 = kf.Phikk_1 * kf.Pk * kf.Phikk_1';
%              kf.Zk = IMU_ven(i,:)';
%              kf.Kk = kf.Pkk_1*kf.Hk'*invbc(kf.Hk*kf.Pkk_1*kf.Hk'+kf.Rk);
% %                 kf.Kk = kf.Pkk_1*kf.Hk'*invbc(kf.Hk*kf.Pkk_1*kf.Hk');
%              kf.Xk = kf.Xkk_1 + kf.Kk*(kf.Zk - kf.Hk*kf.Xkk_1);
%              kf.Pk = (eye(kf.n) - kf.Kk*kf.Hk)*kf.Pkk_1*(eye(kf.n) - kf.Kk*kf.Hk)' + kf.Kk*kf.Rk*kf.Kk';
% %  kf.Pk = (eye(kf.n) - kf.Kk*kf.Hk)*kf.Pkk_1*(eye(kf.n) - kf.Kk*kf.Hk)' ;
%              
% %              kf.Xk = kf.Phikk_1 * kf.Xk;
% %              kf.Pk = kf.Phikk_1 * kf.Pk * kf.Phikk_1' + kf.Qk;
% 
%              %在Kalman更新完成后，利用获取的状态误差，更新 IMU解算值       
%            %%
%              %防止过度收敛是限制方差阵，反馈是为了保持线性度，部分反馈是防止剧烈振荡
%              %滤波状态保留了90％，对状态修正了10%，采纳是全部采纳，因为P阵并没有改变
% 
% %              qnb = qdelphi(qnb, 0.1*kf.Xk(4:6,1));
% %              kf.Xk(4:6,1) = 0.9*kf.Xk(4:6,1);
% %              IMU_ven(i,:) = IMU_ven(i,:) - (0.1*kf.Xk(1:3,1))';
% %              kf.Xk(1:3,1) = 0.9*kf.Xk(1:3,1);  
% 
%              qnb = qdelphi(qnb, kf.Xk(4:6,1));
%              kf.Xk(4:6,1) = kf.Xk(4:6,1);
%              IMU_ven(i,:) = IMU_ven(i,:) - (kf.Xk(1:3,1))';
%              kf.Xk(1:3,1) = kf.Xk(1:3,1);  
% 
%              IMU_att(i,:) = q2att(qnb)'; 
%              Cnb = q2mat(qnb);
% %              xkpk(fix(i/nn),:) = [kf.Xk; diag(kf.Pk)]';   
%              xkpk(fix(i/nn),:) = [kf.Xk; diag(kf.Pk)]'; 
%              
%              IMU_kfattk(mm,:) = IMU_att(i,:);
%              IMU_kfvnk(mm,:) = IMU_ven(i,:);
%              mm = mm+1;
%              
%          
%              
%              
%         end                
%        %% 正常IMU惯导解算 
%         %更新姿态         
%         delphi = imu(i,1:3)'*ts;               %姿态增量
%         qnb = qdelphi(qnb, delphi);                 %更新i时刻的四元数
%         IMU_att(i,:) = q2att(qnb)';                  %获取i时刻的姿态 弧度
%         Cnb = q2mat(qnb);                           %更新i时刻的姿态矩阵    
%         %更新速度 位置        
%         dven = (Cnb*imu(i,4:6)'/ts+eth.gn)*ts;    %获取i时刻的加速度
%         if i < len
%             IMU_ven(i+1,:) = IMU_ven(i,:) + dven';       %获取i+1时刻的速度
%         end
%     end
% %      IMU_att(:,len)/glv.deg
% %      IMU_att = IMU_att/glv.deg;
%     
% 
% 
% 
% 
% 
% disp('自己_kalman对准结果：');   
% IMU_att(len,:)/glv.deg
% 
% [att0v, attkv,vnk, xkpksb] = alignvn(imu, attsb, avp0(7:9), phi0, imuerr, wvn);
% disp('YGM_kalman对准结果：');   
% attkv(length(attkv),:)/glv.deg
% 
% [att0v1, attkv1,vnk1, xkpksb1] = alignvn1(imu, attsb, avp0(7:9), phi0, imuerr, wvn);
% disp('YGM_加补偿kalman对准结果：');   
% attkv1(length(attkv1),:)/glv.deg
% 
% 
% 
% %% 滤波结果的比较
%     t = (1:fix(len/nn))'*nts;
%     %对准姿态
% %     myfigure,
% %     subplot(311);  xygo('phiE/度')
% %     plot(t,[IMU_kfattk(:,1),attkv(:,1),attkv1(:,1)]/glv.deg-3);
% %     legend('My', 'YGM','YGM_compesation');        title('对准姿态');
% %     subplot(312);  xygo('phiN/度');
% %     plot(t,[IMU_kfattk(:,2),attkv(:,2),attkv1(:,2)]/glv.deg-2);
% %     legend('My', 'YGM','YGM_compesation');
% %     subplot(313);  xygo('phiU/度');
% %     plot(t,[IMU_kfattk(:,3),attkv(:,3),attkv1(:,3)]/glv.deg-30);
% %     legend('My', 'YGM','YGM_compesation');
%     
%     %对准速度
% %     myfigure,
% %     subplot(311);  xygo('速度E')
% %     plot(t,[IMU_kfvnk(:,1),vnk(:,1),vnk1(:,1)]);
% %     legend('My', 'YGM','YGM_compesation');    title('对准速度');
% %     subplot(312);  xygo('速度N');
% %     plot(t,[IMU_kfvnk(:,2),vnk(:,2),vnk1(:,2)]);
% %     legend('My', 'YGM','YGM_compesation');
% %     subplot(313);  xygo('速度U');
% %     plot(t,[IMU_kfvnk(:,3),vnk1(:,3),vnk1(:,3)]);
% %     legend('My', 'YGM','YGM_compesation');
%     
%     
%  %自己的Kalman对准姿态误差估计---------------------------------------- ---       
%     myfigure,
%     subplot(331);  xygo('phiE/')
%     plot(t,IMU_kfattk(:,1)/glv.deg);    title('自己的Kalman姿态状态估计及方差');
%     subplot(332);  
%     plot(t,xkpk(:,4)/glv.deg);
%     subplot(333);  
%     plot(t,xkpk(:,16));
%     
%     subplot(334);  xygo('phiN/')
%     plot(t,IMU_kfattk(:,2)/glv.deg);    
%     subplot(335);  
%     plot(t,xkpk(:,5)/glv.deg);
%     subplot(336);  
%     plot(t,xkpk(:,17));
%     
%     subplot(337);  xygo('phiU/')
%     plot(t,IMU_kfattk(:,3)/glv.deg);    
%     subplot(338);  
%     plot(t,xkpk(:,6)/glv.deg);
%     subplot(339);  
%     plot(t,xkpk(:,18));
%     
%     %YGM没有补偿的Kalman对准姿态误差估计--------------------------------------------------        
%     myfigure,
%     subplot(331);  xygo('phiE/')
%     plot(t,attkv(:,1)/glv.deg);    title('YGM的没有补偿的Kalman姿态状态估计及方差');
%     subplot(332);  
%     plot(t,xkpksb(:,4)/glv.deg);
%     subplot(333);  
%     plot(t,xkpksb(:,16));
%     
%     subplot(334);  xygo('phiN/')
%     plot(t,attkv(:,2)/glv.deg);    
%     subplot(335);  
%     plot(t,xkpksb(:,5)/glv.deg);
%     subplot(336);  
%     plot(t,xkpksb(:,17));
%     
%     subplot(337);  xygo('phiU/')
%     plot(t,attkv(:,3)/glv.deg);    
%     subplot(338);  
%     plot(t,xkpksb(:,6)/glv.deg);
%     subplot(339);  
%     plot(t,xkpksb(:,18));
%     
% %YGM有补偿的Kalman对准姿态误差估计--------------------------------------------------        
%     myfigure,
%     subplot(331);  xygo('phiE/')
%     plot(t,attkv1(:,1)/glv.deg);    title('YGM的有补偿的Kalman姿态状态估计及方差');
%     subplot(332);  
%     plot(t,xkpksb1(:,4)/glv.deg);
%     subplot(333);  
%     plot(t,xkpksb1(:,16));
%     
%     subplot(334);  xygo('phiN/')
%     plot(t,attkv1(:,2)/glv.deg);    
%     subplot(335);  
%     plot(t,xkpksb1(:,5)/glv.deg);
%     subplot(336);  
%     plot(t,xkpksb1(:,17));
%     
%     subplot(337);  xygo('phiU/')
%     plot(t,attkv1(:,3)/glv.deg);    
%     subplot(338);  
%     plot(t,xkpksb1(:,6)/glv.deg);
%     subplot(339);  
%     plot(t,xkpksb1(:,18));
% % 
% % 
% %     
% %  %自己的Kalman对准的  速度 误差估计---------------------------------------- ---       
% %     myfigure,
% %     subplot(331);  xygo('E/')
% %     plot(t,IMU_kfvnk(:,1));    title('自己的Kalman计算速度值');
% %     subplot(332);  
% %     plot(t,xkpk(:,1));   title('自己的Kalman计算误差值');
% %     subplot(333);  
% %     plot(t,xkpk(:,13));   title('自己的Kalman计算方差值');
% %     
% %     subplot(334);  xygo('N/')
% %     plot(t,IMU_kfvnk(:,2));    
% %     subplot(335);  
% %     plot(t,xkpk(:,2));
% %     subplot(336);  
% %     plot(t,xkpk(:,14));
% %     
% %     subplot(337);  xygo('U/')
% %     plot(t,IMU_kfvnk(:,3));    
% %     subplot(338);  
% %     plot(t,xkpk(:,3));
% %     subplot(339);  
% %     plot(t,xkpk(:,15));
% %     
% % %YGM有补偿的Kalman对准的 速度 误差估计--------------------------------------------------        
% %     myfigure,
% %     subplot(331);  xygo('E/')
% %     plot(t,vnk1(:,1));    title('YGM的有补偿的Kalman计算速度值');
% %     subplot(332);  
% %     plot(t,xkpksb1(:,1));  title('YGM的有补偿的Kalman计算误差值');
% %     subplot(333);  
% %     plot(t,xkpksb1(:,13));  title('YGM的有补偿的Kalman计算方差值');
% %     
% %     subplot(334);  xygo('N/')
% %     plot(t,vnk1(:,2));    
% %     subplot(335);  
% %     plot(t,xkpksb1(:,2));
% %     subplot(336);  
% %     plot(t,xkpksb1(:,14));
% %     
% %     subplot(337);  xygo('U/')
% %     plot(t,vnk1(:,3));    
% %     subplot(338);  
% %     plot(t,xkpksb1(:,3));
% %     subplot(339);  
% %     plot(t,xkpksb1(:,15));
% % %YGM没有补偿的Kalman对准的 速度 误差估计--------------------------------------------------        
%     myfigure,
%     subplot(331);  xygo('E/')
%     plot(t,vnk(:,1));    title('YGM的没有补偿的Kalman计算速度值');
%     subplot(332);  
%     plot(t,xkpksb(:,1));  title('YGM的没有补偿的Kalman计算误差值');
%     subplot(333);  
%     plot(t,xkpksb(:,13));  title('YGM的没有补偿的Kalman计算方差值');
%     
%     subplot(334);  xygo('N/')
%     plot(t,vnk(:,2));    
%     subplot(335);  
%     plot(t,xkpksb(:,2));
%     subplot(336);  
%     plot(t,xkpksb(:,14));
%     
%     subplot(337);  xygo('U/')
%     plot(t,vnk(:,3));    
%     subplot(338);  
%     plot(t,xkpksb(:,3));
%     subplot(339);  
%     plot(t,xkpksb(:,15));