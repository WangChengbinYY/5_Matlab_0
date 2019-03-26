function kf = kfinit_wcb(n,m)
% ��������׼���ٶȷ���Kalman���� ��ʼ��
% 2017.11.27 ���ɱ�

kf = [];
kf.n = n; % ״̬����ά��
kf.m = m; % �۲ⷽ��ά��
kf.I = eye(n);  %��λ����
kf.Xk_1 = zeros(n,1);
kf.Pkk_1 = zeros(n,n);
kf.Xk = zeros(n,1);
kf.Pk = zeros(n,n);
kf.Zk = zeros(m,1);
kf.Hk = zeros(m,n);
kf.Kk = zeros(n,m);
kf.XkZk = zeros(n,1);  %��¼�м�� ����Zk��Xk�Ĺ��� Kk(Zk-HkXkk_1);

kf.Ft = zeros(n,n);     % ����ϵͳ��״̬����
kf.Qt = zeros(kf.n);    % ����ϵͳ��״̬����ϵͳ����������
kf.Phikk_1 = zeros(n,n);   % ��ɢ����״̬����
kf.Qk = zeros(n,n);
kf.Rk = zeros(m,m);



