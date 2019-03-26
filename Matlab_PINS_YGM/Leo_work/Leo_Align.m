% 10/20/2017 Modified by Chengbin Wang, Tsinghua University
%
%
%% 
%--------------------------------------------
% 初始化全局变量 glv，并获取仿真数据 
% imu 包括 陀螺 加计 时间 数据
% parameter 包括 各种仿真参数
%--------------------------------------------
    clear all;
    glvs
    [rpath, dpath, mytestflag] = psinsenvi();
    dpath = [dpath, '\TestData_Align.mat'];
    load(dpath);
    disp(describe);
insupdate

%%  
%------------------------------------------
% 先利用1分钟的陀螺和加计数据进行 静基座 粗对准 获取初始姿态 attsb
%------------------------------------------
    num = fix(60/ts);
    imu_sb = imu(1:num,:);
    [attsb, qnb] = alignsb(imu, avp0(7:9));
    disp('静基座1分钟粗对准结果：');   
    attsb/glv.deg

%%
%------------------------------------------
% 在静基座粗对准基础上开始 速度反馈的Kalman滤波精对准（不包含圆锥 划船补偿）
%------------------------------------------

%% ----惯导解算初始化-----

    nn = 4; nts = nn*ts;                        %没隔nn个采样周期进行一次Kalman滤波
    len = fix(length(imu)/nn)*nn;          %数据长度为len
    [IMU_att,IMU_ven] = prealloc(len, 3, 3);    %惯导解算出来的姿态和速度数据 姿态单位为弧度
    phi0 = [0.5;0.5;2]*glv.deg;                     %设定初始失准角
    IMU_att(1,:) = attsb'; 
    wvn = [0.01;0.01;0.01];                  %粗 设定初始速度    
    IMU_ven(1,:) = wvn';
    Cnb = a2mat(attsb); qnb = a2qua(attsb);

%% ----Kalman滤波初始化-----
    eth = earth(avp0(7:9));                 %获取地球相关计算参数
    v_error = zeros(3,1);  p_error = zeros(3,1); %速度、位置误差
%第一步 进行陀螺、加计参数设定
    % imuerr;
%第二步 进行Kalmn滤波参数初始设定 包括Ft0,Qt0, Rt0, Pt0, Ht0及初始姿态
    kf = [];
    kf.nts = nts;
    kf.Cnb = a2mat(attsb(1:3));     %利用粗对准结果
    kf.qnb = a2qua(attsb(1:3));
    kf.n = 12; kf.m = 3;            %n 状态量个数，m 观测量个数
    kf.I = eye(kf.n);
    kf.Xk = zeros(kf.n,1); kf.Zk = zeros(kf.m,1);
    kf.Ft = zeros(kf.n); kf.Hk = zeros(kf.m,kf.n);
    kf.Qk = zeros(kf.n); kf.Qt = zeros(kf.n); 
    kf.Rk = zeros(kf.m);  kf.Pk = zeros(kf.n); 
    kf.Pkk_1 = zeros(kf.n); kf.Xkk_1 = zeros(kf.n,1); 
    kf.Kk = zeros(kf.n);  kf.Phikk_1 = zeros(kf.n);

%     kf.Xk = [wvn; phi0; [0;0;0]; [0;0;0]];   %系统初始状态 全设为0
    kf.Pk = diag([[1;1;1]; phi0; imuerr.eb; imuerr.db])^2;      %初始 P0
    kf.Rk = diag(IMU_ven(1,1))^2;

    [xkpk,xkpksb,xkpksb1] = prealloc(fix(len/nn), 2*kf.n, 2*kf.n, 2*kf.n);  %将中间系统状态 数据记录下来
    [IMU_kfattk,IMU_kfvnk] = prealloc(fix(len/nn), 3, 3); 

    %先计算一下 Ft中的不变量
    [wie,lat,RMh,RNh] = ... % just for short
    setvals(glv.wie,avp0(7,1),eth.RMh+avp0(9,1),eth.RNh+avp0(9,1)); 
    kf.Ft(1,2) = 2*wie*sin(lat);        kf.Ft(2,1) = -kf.Ft(1,2);
    kf.Ft(1,3) = -2*wie*cos(lat);       kf.Ft(3,1) = -kf.Ft(1,3);   
    kf.Ft(1,5) = -eth.g;                kf.Ft(2,4) = eth.g;
    kf.Ft(4,2) = -1/RMh;                kf.Ft(5,1) = 1/RNh;         kf.Ft(6,1) = tan(lat)/RNh; %后面可以试验，看影响大不！！？？
    kf.Ft(4,5) = wie*sin(lat);          kf.Ft(5,4) = -kf.Ft(4,5);
    kf.Ft(4,6) = wie*cos(lat);          kf.Ft(6,4) = -kf.Ft(4,6);
    kf.Ft(1:3,10:12) = kf.Cnb;          kf.Ft(4:6,7:9) = -kf.Cnb;      %是否取Cnb 或者 eye(3) 计算没有影响！！！！！！
%     kf.Ft(1:3,10:12) = eye(3);          kf.Ft(4:6,7:9) = -eye(3);
    kf.Hk(:,1:3) = eye(3);
    
     [kfZk_HkXkk_1,IMU_Xk_add] = prealloc(fix(len/nn), 6, 6); 
%% ----开始循环惯导解算，每隔nn个采样周期 加入一次Kalman滤波-----
    dven = zeros(3,1);  delphi = zeros(3,1);  mm = 1;
    for i=1:1:len
        %% 判断若是满足 nn 间隔，则加入 Kalman滤波  
        if mod(i,nn) == 0
             %计算新时刻的Ft
              kf.Ft(1:3,10:12) = kf.Cnb; kf.Ft(4:6,7:9) = -kf.Cnb;
             %计算新时刻的Qt
             % 按照陀螺 加计的随机游走系数web wdb 来计算 方差强度阵q
             % ！！！此处的web wdb 随机游走系数已经做过 单位换算了！！！,并且直接和系统状态方程相关，不通用
             kf.Qt = Leo_GetQtFrom(kf.n,Cnb,imuerr.web,imuerr.wdb);  %后面做实验 直接用固定的，是否有变化             
             [kf.Phikk_1, kf.Qk] = kfc2d(kf.Ft, kf.Qt, ts*nn, 3);  
            
             kf.Xkk_1 = kf.Phikk_1 * kf.Xk;
             kf.Pkk_1 = kf.Phikk_1 * kf.Pk * kf.Phikk_1' + kf.Qk;
%                 kf.Pkk_1 = kf.Phikk_1 * kf.Pk * kf.Phikk_1';
             kf.Zk = IMU_ven(i,:)';
             kf.Kk = kf.Pkk_1*kf.Hk'*invbc(kf.Hk*kf.Pkk_1*kf.Hk'+kf.Rk);
%                 kf.Kk = kf.Pkk_1*kf.Hk'*invbc(kf.Hk*kf.Pkk_1*kf.Hk');
             kf.Xk = kf.Xkk_1 + kf.Kk*(kf.Zk - kf.Hk*kf.Xkk_1);
             kf.Pk = (eye(kf.n) - kf.Kk*kf.Hk)*kf.Pkk_1*(eye(kf.n) - kf.Kk*kf.Hk)' + kf.Kk*kf.Rk*kf.Kk';
%  kf.Pk = (eye(kf.n) - kf.Kk*kf.Hk)*kf.Pkk_1*(eye(kf.n) - kf.Kk*kf.Hk)' ;
             
%              kf.Xk = kf.Phikk_1 * kf.Xk;
%              kf.Pk = kf.Phikk_1 * kf.Pk * kf.Phikk_1' + kf.Qk;

             %在Kalman更新完成后，利用获取的状态误差，更新 IMU解算值       
           %%
             %防止过度收敛是限制方差阵，反馈是为了保持线性度，部分反馈是防止剧烈振荡
             %滤波状态保留了90％，对状态修正了10%，采纳是全部采纳，因为P阵并没有改变

%              qnb = qdelphi(qnb, 0.1*kf.Xk(4:6,1));
%              kf.Xk(4:6,1) = 0.9*kf.Xk(4:6,1);
%              IMU_ven(i,:) = IMU_ven(i,:) - (0.1*kf.Xk(1:3,1))';
%              kf.Xk(1:3,1) = 0.9*kf.Xk(1:3,1);  

             qnb = qdelphi(qnb, kf.Xk(4:6,1));
             kf.Xk(4:6,1) = kf.Xk(4:6,1);
             IMU_ven(i,:) = IMU_ven(i,:) - (kf.Xk(1:3,1))';
             kf.Xk(1:3,1) = kf.Xk(1:3,1);  

             IMU_att(i,:) = q2att(qnb)'; 
             Cnb = q2mat(qnb);
%              xkpk(fix(i/nn),:) = [kf.Xk; diag(kf.Pk)]';   
             xkpk(fix(i/nn),:) = [kf.Xk; diag(kf.Pk)]'; 
             
             IMU_kfattk(mm,:) = IMU_att(i,:);
             IMU_kfvnk(mm,:) = IMU_ven(i,:);
             mm = mm+1;
             
         
             
             
        end                
       %% 正常IMU惯导解算 
        %更新姿态         
        delphi = imu(i,1:3)'*ts;               %姿态增量
        qnb = qdelphi(qnb, delphi);                 %更新i时刻的四元数
        IMU_att(i,:) = q2att(qnb)';                  %获取i时刻的姿态 弧度
        Cnb = q2mat(qnb);                           %更新i时刻的姿态矩阵    
        %更新速度 位置        
        dven = (Cnb*imu(i,4:6)'/ts+eth.gn)*ts;    %获取i时刻的加速度
        if i < len
            IMU_ven(i+1,:) = IMU_ven(i,:) + dven';       %获取i+1时刻的速度
        end
    end
%      IMU_att(:,len)/glv.deg
%      IMU_att = IMU_att/glv.deg;
    





disp('自己_kalman对准结果：');   
IMU_att(len,:)/glv.deg

[att0v, attkv,vnk, xkpksb] = alignvn(imu, attsb, avp0(7:9), phi0, imuerr, wvn);
disp('YGM_kalman对准结果：');   
attkv(length(attkv),:)/glv.deg

[att0v1, attkv1,vnk1, xkpksb1] = alignvn1(imu, attsb, avp0(7:9), phi0, imuerr, wvn);
disp('YGM_加补偿kalman对准结果：');   
attkv1(length(attkv1),:)/glv.deg



%% 滤波结果的比较
    t = (1:fix(len/nn))'*nts;
    %对准姿态
%     myfigure,
%     subplot(311);  xygo('phiE/度')
%     plot(t,[IMU_kfattk(:,1),attkv(:,1),attkv1(:,1)]/glv.deg-3);
%     legend('My', 'YGM','YGM_compesation');        title('对准姿态');
%     subplot(312);  xygo('phiN/度');
%     plot(t,[IMU_kfattk(:,2),attkv(:,2),attkv1(:,2)]/glv.deg-2);
%     legend('My', 'YGM','YGM_compesation');
%     subplot(313);  xygo('phiU/度');
%     plot(t,[IMU_kfattk(:,3),attkv(:,3),attkv1(:,3)]/glv.deg-30);
%     legend('My', 'YGM','YGM_compesation');
    
    %对准速度
%     myfigure,
%     subplot(311);  xygo('速度E')
%     plot(t,[IMU_kfvnk(:,1),vnk(:,1),vnk1(:,1)]);
%     legend('My', 'YGM','YGM_compesation');    title('对准速度');
%     subplot(312);  xygo('速度N');
%     plot(t,[IMU_kfvnk(:,2),vnk(:,2),vnk1(:,2)]);
%     legend('My', 'YGM','YGM_compesation');
%     subplot(313);  xygo('速度U');
%     plot(t,[IMU_kfvnk(:,3),vnk1(:,3),vnk1(:,3)]);
%     legend('My', 'YGM','YGM_compesation');
    
    
 %自己的Kalman对准姿态误差估计---------------------------------------- ---       
    myfigure,
    subplot(331);  xygo('phiE/')
    plot(t,IMU_kfattk(:,1)/glv.deg);    title('自己的Kalman姿态状态估计及方差');
    subplot(332);  
    plot(t,xkpk(:,4)/glv.deg);
    subplot(333);  
    plot(t,xkpk(:,16));
    
    subplot(334);  xygo('phiN/')
    plot(t,IMU_kfattk(:,2)/glv.deg);    
    subplot(335);  
    plot(t,xkpk(:,5)/glv.deg);
    subplot(336);  
    plot(t,xkpk(:,17));
    
    subplot(337);  xygo('phiU/')
    plot(t,IMU_kfattk(:,3)/glv.deg);    
    subplot(338);  
    plot(t,xkpk(:,6)/glv.deg);
    subplot(339);  
    plot(t,xkpk(:,18));
    
    %YGM没有补偿的Kalman对准姿态误差估计--------------------------------------------------        
    myfigure,
    subplot(331);  xygo('phiE/')
    plot(t,attkv(:,1)/glv.deg);    title('YGM的没有补偿的Kalman姿态状态估计及方差');
    subplot(332);  
    plot(t,xkpksb(:,4)/glv.deg);
    subplot(333);  
    plot(t,xkpksb(:,16));
    
    subplot(334);  xygo('phiN/')
    plot(t,attkv(:,2)/glv.deg);    
    subplot(335);  
    plot(t,xkpksb(:,5)/glv.deg);
    subplot(336);  
    plot(t,xkpksb(:,17));
    
    subplot(337);  xygo('phiU/')
    plot(t,attkv(:,3)/glv.deg);    
    subplot(338);  
    plot(t,xkpksb(:,6)/glv.deg);
    subplot(339);  
    plot(t,xkpksb(:,18));
    
%YGM有补偿的Kalman对准姿态误差估计--------------------------------------------------        
    myfigure,
    subplot(331);  xygo('phiE/')
    plot(t,attkv1(:,1)/glv.deg);    title('YGM的有补偿的Kalman姿态状态估计及方差');
    subplot(332);  
    plot(t,xkpksb1(:,4)/glv.deg);
    subplot(333);  
    plot(t,xkpksb1(:,16));
    
    subplot(334);  xygo('phiN/')
    plot(t,attkv1(:,2)/glv.deg);    
    subplot(335);  
    plot(t,xkpksb1(:,5)/glv.deg);
    subplot(336);  
    plot(t,xkpksb1(:,17));
    
    subplot(337);  xygo('phiU/')
    plot(t,attkv1(:,3)/glv.deg);    
    subplot(338);  
    plot(t,xkpksb1(:,6)/glv.deg);
    subplot(339);  
    plot(t,xkpksb1(:,18));
% 
% 
%     
%  %自己的Kalman对准的  速度 误差估计---------------------------------------- ---       
%     myfigure,
%     subplot(331);  xygo('E/')
%     plot(t,IMU_kfvnk(:,1));    title('自己的Kalman计算速度值');
%     subplot(332);  
%     plot(t,xkpk(:,1));   title('自己的Kalman计算误差值');
%     subplot(333);  
%     plot(t,xkpk(:,13));   title('自己的Kalman计算方差值');
%     
%     subplot(334);  xygo('N/')
%     plot(t,IMU_kfvnk(:,2));    
%     subplot(335);  
%     plot(t,xkpk(:,2));
%     subplot(336);  
%     plot(t,xkpk(:,14));
%     
%     subplot(337);  xygo('U/')
%     plot(t,IMU_kfvnk(:,3));    
%     subplot(338);  
%     plot(t,xkpk(:,3));
%     subplot(339);  
%     plot(t,xkpk(:,15));
%     
% %YGM有补偿的Kalman对准的 速度 误差估计--------------------------------------------------        
%     myfigure,
%     subplot(331);  xygo('E/')
%     plot(t,vnk1(:,1));    title('YGM的有补偿的Kalman计算速度值');
%     subplot(332);  
%     plot(t,xkpksb1(:,1));  title('YGM的有补偿的Kalman计算误差值');
%     subplot(333);  
%     plot(t,xkpksb1(:,13));  title('YGM的有补偿的Kalman计算方差值');
%     
%     subplot(334);  xygo('N/')
%     plot(t,vnk1(:,2));    
%     subplot(335);  
%     plot(t,xkpksb1(:,2));
%     subplot(336);  
%     plot(t,xkpksb1(:,14));
%     
%     subplot(337);  xygo('U/')
%     plot(t,vnk1(:,3));    
%     subplot(338);  
%     plot(t,xkpksb1(:,3));
%     subplot(339);  
%     plot(t,xkpksb1(:,15));
% %YGM没有补偿的Kalman对准的 速度 误差估计--------------------------------------------------        
    myfigure,
    subplot(331);  xygo('E/')
    plot(t,vnk(:,1));    title('YGM的没有补偿的Kalman计算速度值');
    subplot(332);  
    plot(t,xkpksb(:,1));  title('YGM的没有补偿的Kalman计算误差值');
    subplot(333);  
    plot(t,xkpksb(:,13));  title('YGM的没有补偿的Kalman计算方差值');
    
    subplot(334);  xygo('N/')
    plot(t,vnk(:,2));    
    subplot(335);  
    plot(t,xkpksb(:,2));
    subplot(336);  
    plot(t,xkpksb(:,14));
    
    subplot(337);  xygo('U/')
    plot(t,vnk(:,3));    
    subplot(338);  
    plot(t,xkpksb(:,3));
    subplot(339);  
    plot(t,xkpksb(:,15));