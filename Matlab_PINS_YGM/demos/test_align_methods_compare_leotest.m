% Several popular align methods are compared.
% See also  test_align_methods_compare_lgimu, test_align_ekf, test_align_ukf,
%           alignvn, aligni0, aligncmps.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 10/03/2014
glvs
ts = 0.1;   % sampling interval
T = 1000;
avp0 = avpset([0;0;30], [0;0;0], [30;108;380]);
imuerr = imuerrset(0.01, 50, 0.0001, 0.10);
imu = imustatic(avp0, ts, T, imuerr);   % IMU simulation
davp = avpseterr([-30;30;30], [0.01;0.01;0.01]*0, [1;1;1]*0);
avp = avpadderr(avp0, davp);

%% coarse align
[attsb, qnb] = alignsb(imu, avp(7:9));
attsb/glv.deg
%% vn-meas. Kalman filter
phi = [.5;.5;5]*glv.deg;
wvn = [0.01;0.01;0.01];
[att0v, attkv, xkpkv] = alignvn1(imu, qnb, avp(7:9), phi, imuerr, wvn);
att0v/glv.deg
%% gyro-compass method
ctl0 = [20; 30]; ctl1 = [50; 300];
[att0c, attkc] = aligncmps(imu, qnb, avp0(7:9), ctl0, ctl1);
att0c/glv.deg

t = (1:length(attkv))'*ts*4;
plot(t, [attkv(:,3),attkc(:,3)]/glv.deg)
legend( 'Kalman vn', 'gyro-compass');
%% compare & show different methods
% [phii0p,phii0v,phii0w,phikf,phikv,phic] = setvals(zeros(size(attkc)));
% for k=1:length(attkc)
%     phii0p(k,:) = aa2phi(resi0.attk(k,:)',avp0(1:3))';
%     phii0v(k,:) = aa2phi(resi0.attkv(k,:)',avp0(1:3))';
%     phii0w(k,:) = aa2phi(attkw(k,:)',avp0(1:3))';
%     phikf(k,:) = aa2phi(attkf(k,:)',avp0(1:3))';
%     phikv(k,:) = aa2phi(attkv(k,:)',avp0(1:3))';
%     phic(k,:) = aa2phi(attkc(k,:)',avp0(1:3))';
% end
% t = (1:length(phic))'*resi0.nts;
% myfigure,
% subplot(211);  xygo('phiE')
% plot(t, [phii0p(:,1),phii0v(:,1),phii0w(:,1),phikf(:,1),phikv(:,1),phic(:,1)]/glv.sec)
% legend('i0 pos', 'i0 vel', 'i0 Whaba', 'Kalman fn', 'Kalman vn', 'gyro-compass');
% subplot(212);  xygo('phiN')
% plot(t, [phii0p(:,2),phii0v(:,2),phii0w(:,2),phikf(:,2),phikv(:,2),phic(:,2)]/glv.sec)
% legend('i0 pos', 'i0 vel', 'i0 Whaba', 'Kalman fn', 'Kalman vn', 'gyro-compass');
% myfigure, xygo('phiU')
% plot(t, [phii0p(:,3),phii0v(:,3),phii0w(:,3),phikf(:,3),phikv(:,3),phic(:,3)]/glv.min)
% legend('i0 pos', 'i0 vel', 'i0 Whaba', 'Kalman fn', 'Kalman vn', 'gyro-compass');