% 10/20/2017 Modified by Chengbin Wang, Tsinghua University
%
%
%% 
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
    disp(describe);
insupdate

%%  
%------------------------------------------
% ������1���ӵ����ݺͼӼ����ݽ��� ������ �ֶ�׼ ��ȡ��ʼ��̬ attsb
%------------------------------------------
    num = fix(60/ts);
    imu_sb = imu(1:num,:);
    [attsb, qnb] = alignsb(imu, avp0(7:9));
    disp('������1���Ӵֶ�׼�����');   
    attsb/glv.deg

%%
%------------------------------------------
% �ھ������ֶ�׼�����Ͽ�ʼ �ٶȷ�����Kalman�˲�����׼��������Բ׶ ����������
%------------------------------------------

%% ----�ߵ������ʼ��-----

    nn = 4; nts = nn*ts;                        %û��nn���������ڽ���һ��Kalman�˲�
    len = fix(length(imu)/nn)*nn;          %���ݳ���Ϊlen
    [IMU_att,IMU_ven] = prealloc(len, 3, 3);    %�ߵ������������̬���ٶ����� ��̬��λΪ����
    phi0 = [0.5;0.5;2]*glv.deg;                     %�趨��ʼʧ׼��
    IMU_att(1,:) = attsb'; 
    wvn = [0.01;0.01;0.01];                  %�� �趨��ʼ�ٶ�    
    IMU_ven(1,:) = wvn';
    Cnb = a2mat(attsb); qnb = a2qua(attsb);

%% ----Kalman�˲���ʼ��-----
    eth = earth(avp0(7:9));                 %��ȡ������ؼ������
    v_error = zeros(3,1);  p_error = zeros(3,1); %�ٶȡ�λ�����
%��һ�� �������ݡ��ӼƲ����趨
    % imuerr;
%�ڶ��� ����Kalmn�˲�������ʼ�趨 ����Ft0,Qt0, Rt0, Pt0, Ht0����ʼ��̬
    kf = [];
    kf.nts = nts;
    kf.Cnb = a2mat(attsb(1:3));     %���ôֶ�׼���
    kf.qnb = a2qua(attsb(1:3));
    kf.n = 12; kf.m = 3;            %n ״̬��������m �۲�������
    kf.I = eye(kf.n);
    kf.Xk = zeros(kf.n,1); kf.Zk = zeros(kf.m,1);
    kf.Ft = zeros(kf.n); kf.Hk = zeros(kf.m,kf.n);
    kf.Qk = zeros(kf.n); kf.Qt = zeros(kf.n); 
    kf.Rk = zeros(kf.m);  kf.Pk = zeros(kf.n); 
    kf.Pkk_1 = zeros(kf.n); kf.Xkk_1 = zeros(kf.n,1); 
    kf.Kk = zeros(kf.n);  kf.Phikk_1 = zeros(kf.n);

%     kf.Xk = [wvn; phi0; [0;0;0]; [0;0;0]];   %ϵͳ��ʼ״̬ ȫ��Ϊ0
    kf.Pk = diag([[1;1;1]; phi0; imuerr.eb; imuerr.db])^2;      %��ʼ P0
    kf.Rk = diag(IMU_ven(1,1))^2;

    [xkpk,xkpksb,xkpksb1] = prealloc(fix(len/nn), 2*kf.n, 2*kf.n, 2*kf.n);  %���м�ϵͳ״̬ ���ݼ�¼����
    [IMU_kfattk,IMU_kfvnk] = prealloc(fix(len/nn), 3, 3); 

    %�ȼ���һ�� Ft�еĲ�����
    [wie,lat,RMh,RNh] = ... % just for short
    setvals(glv.wie,avp0(7,1),eth.RMh+avp0(9,1),eth.RNh+avp0(9,1)); 
    kf.Ft(1,2) = 2*wie*sin(lat);        kf.Ft(2,1) = -kf.Ft(1,2);
    kf.Ft(1,3) = -2*wie*cos(lat);       kf.Ft(3,1) = -kf.Ft(1,3);   
    kf.Ft(1,5) = -eth.g;                kf.Ft(2,4) = eth.g;
    kf.Ft(4,2) = -1/RMh;                kf.Ft(5,1) = 1/RNh;         kf.Ft(6,1) = tan(lat)/RNh; %����������飬��Ӱ��󲻣�������
    kf.Ft(4,5) = wie*sin(lat);          kf.Ft(5,4) = -kf.Ft(4,5);
    kf.Ft(4,6) = wie*cos(lat);          kf.Ft(6,4) = -kf.Ft(4,6);
    kf.Ft(1:3,10:12) = kf.Cnb;          kf.Ft(4:6,7:9) = -kf.Cnb;      %�Ƿ�ȡCnb ���� eye(3) ����û��Ӱ�죡����������
%     kf.Ft(1:3,10:12) = eye(3);          kf.Ft(4:6,7:9) = -eye(3);
    kf.Hk(:,1:3) = eye(3);
    
     [kfZk_HkXkk_1,IMU_Xk_add] = prealloc(fix(len/nn), 6, 6); 
%% ----��ʼѭ���ߵ����㣬ÿ��nn���������� ����һ��Kalman�˲�-----
    dven = zeros(3,1);  delphi = zeros(3,1);  mm = 1;
    for i=1:1:len
        %% �ж��������� nn ���������� Kalman�˲�  
        if mod(i,nn) == 0
             %������ʱ�̵�Ft
              kf.Ft(1:3,10:12) = kf.Cnb; kf.Ft(4:6,7:9) = -kf.Cnb;
             %������ʱ�̵�Qt
             % �������� �ӼƵ��������ϵ��web wdb ������ ����ǿ����q
             % �������˴���web wdb �������ϵ���Ѿ����� ��λ�����ˣ�����,����ֱ�Ӻ�ϵͳ״̬������أ���ͨ��
             kf.Qt = Leo_GetQtFrom(kf.n,Cnb,imuerr.web,imuerr.wdb);  %������ʵ�� ֱ���ù̶��ģ��Ƿ��б仯             
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

             %��Kalman������ɺ����û�ȡ��״̬������ IMU����ֵ       
           %%
             %��ֹ�������������Ʒ����󣬷�����Ϊ�˱������Զȣ����ַ����Ƿ�ֹ������
             %�˲�״̬������90������״̬������10%��������ȫ�����ɣ���ΪP��û�иı�

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
       %% ����IMU�ߵ����� 
        %������̬         
        delphi = imu(i,1:3)'*ts;               %��̬����
        qnb = qdelphi(qnb, delphi);                 %����iʱ�̵���Ԫ��
        IMU_att(i,:) = q2att(qnb)';                  %��ȡiʱ�̵���̬ ����
        Cnb = q2mat(qnb);                           %����iʱ�̵���̬����    
        %�����ٶ� λ��        
        dven = (Cnb*imu(i,4:6)'/ts+eth.gn)*ts;    %��ȡiʱ�̵ļ��ٶ�
        if i < len
            IMU_ven(i+1,:) = IMU_ven(i,:) + dven';       %��ȡi+1ʱ�̵��ٶ�
        end
    end
%      IMU_att(:,len)/glv.deg
%      IMU_att = IMU_att/glv.deg;
    





disp('�Լ�_kalman��׼�����');   
IMU_att(len,:)/glv.deg

[att0v, attkv,vnk, xkpksb] = alignvn(imu, attsb, avp0(7:9), phi0, imuerr, wvn);
disp('YGM_kalman��׼�����');   
attkv(length(attkv),:)/glv.deg

[att0v1, attkv1,vnk1, xkpksb1] = alignvn1(imu, attsb, avp0(7:9), phi0, imuerr, wvn);
disp('YGM_�Ӳ���kalman��׼�����');   
attkv1(length(attkv1),:)/glv.deg



%% �˲�����ıȽ�
    t = (1:fix(len/nn))'*nts;
    %��׼��̬
%     myfigure,
%     subplot(311);  xygo('phiE/��')
%     plot(t,[IMU_kfattk(:,1),attkv(:,1),attkv1(:,1)]/glv.deg-3);
%     legend('My', 'YGM','YGM_compesation');        title('��׼��̬');
%     subplot(312);  xygo('phiN/��');
%     plot(t,[IMU_kfattk(:,2),attkv(:,2),attkv1(:,2)]/glv.deg-2);
%     legend('My', 'YGM','YGM_compesation');
%     subplot(313);  xygo('phiU/��');
%     plot(t,[IMU_kfattk(:,3),attkv(:,3),attkv1(:,3)]/glv.deg-30);
%     legend('My', 'YGM','YGM_compesation');
    
    %��׼�ٶ�
%     myfigure,
%     subplot(311);  xygo('�ٶ�E')
%     plot(t,[IMU_kfvnk(:,1),vnk(:,1),vnk1(:,1)]);
%     legend('My', 'YGM','YGM_compesation');    title('��׼�ٶ�');
%     subplot(312);  xygo('�ٶ�N');
%     plot(t,[IMU_kfvnk(:,2),vnk(:,2),vnk1(:,2)]);
%     legend('My', 'YGM','YGM_compesation');
%     subplot(313);  xygo('�ٶ�U');
%     plot(t,[IMU_kfvnk(:,3),vnk1(:,3),vnk1(:,3)]);
%     legend('My', 'YGM','YGM_compesation');
    
    
 %�Լ���Kalman��׼��̬������---------------------------------------- ---       
    myfigure,
    subplot(331);  xygo('phiE/')
    plot(t,IMU_kfattk(:,1)/glv.deg);    title('�Լ���Kalman��̬״̬���Ƽ�����');
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
    
    %YGMû�в�����Kalman��׼��̬������--------------------------------------------------        
    myfigure,
    subplot(331);  xygo('phiE/')
    plot(t,attkv(:,1)/glv.deg);    title('YGM��û�в�����Kalman��̬״̬���Ƽ�����');
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
    
%YGM�в�����Kalman��׼��̬������--------------------------------------------------        
    myfigure,
    subplot(331);  xygo('phiE/')
    plot(t,attkv1(:,1)/glv.deg);    title('YGM���в�����Kalman��̬״̬���Ƽ�����');
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
%  %�Լ���Kalman��׼��  �ٶ� ������---------------------------------------- ---       
%     myfigure,
%     subplot(331);  xygo('E/')
%     plot(t,IMU_kfvnk(:,1));    title('�Լ���Kalman�����ٶ�ֵ');
%     subplot(332);  
%     plot(t,xkpk(:,1));   title('�Լ���Kalman�������ֵ');
%     subplot(333);  
%     plot(t,xkpk(:,13));   title('�Լ���Kalman���㷽��ֵ');
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
% %YGM�в�����Kalman��׼�� �ٶ� ������--------------------------------------------------        
%     myfigure,
%     subplot(331);  xygo('E/')
%     plot(t,vnk1(:,1));    title('YGM���в�����Kalman�����ٶ�ֵ');
%     subplot(332);  
%     plot(t,xkpksb1(:,1));  title('YGM���в�����Kalman�������ֵ');
%     subplot(333);  
%     plot(t,xkpksb1(:,13));  title('YGM���в�����Kalman���㷽��ֵ');
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
% %YGMû�в�����Kalman��׼�� �ٶ� ������--------------------------------------------------        
    myfigure,
    subplot(331);  xygo('E/')
    plot(t,vnk(:,1));    title('YGM��û�в�����Kalman�����ٶ�ֵ');
    subplot(332);  
    plot(t,xkpksb(:,1));  title('YGM��û�в�����Kalman�������ֵ');
    subplot(333);  
    plot(t,xkpksb(:,13));  title('YGM��û�в�����Kalman���㷽��ֵ');
    
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