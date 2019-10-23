function [att0, attk] = alignWahba(imu, pos, ts)
% SINS initial align based on inertial frame and Wahba method.
%
% Prototype: [att0, attk] = alignWahba(imu, pos, ts)
% Inputs: imu - IMU data
%         pos - position
%         ts - IMU sampling interval
% Output: att0 - attitude align result
%
% See also  alignfn, alignvn, aligncmps, aligni0, alignsb.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 29/06/2012
global glv
    if nargin<3,  ts = imu(2,7)-imu(1,7);  end
	nn = 4; nts = nn*ts;
    len = fix(length(imu)/nn)*nn;
	eth = earth(pos);   g0 = -eth.gn(3);
    qib0b = [1; 0; 0; 0];
    [vib0, vi0] = setvals(zeros(3,1));
    K = zeros(4);
    attk = zeros(len/nn, 3);
    ki = timebar(nn, len, 'Initial align using i0 & Wahba method.');
    for k=1:nn:len-nn+1
        wvm = imu(k:k+nn-1, 1:6);  kts = (k+nn-1)*ts;
        [phim, dvbm] = cnscl(wvm);
        vib0 = 0*vib0 + qmulv(qib0b, dvbm);  %加计积分速度
        vi0 = 0*vi0 + [eth.cl*cos(kts*glv.wie); eth.cl*sin(kts*glv.wie); eth.sl]*g0*nts; %重力积分速度
        qib0b = qupdt(qib0b, phim); %旋转矢量更新  b-ib0系
        dM = rq2m([0;vib0])-lq2m([0;vi0]); %论文的式27
        K = 0.99991*K + dM'*dM*nts; %积分函数的更新,初值累加
        [v, d] = eig(K); %返回特征值d和特征向量v
        qi0ib0 = v(:,1); %与一重积分速度求qi0ib0算法不一样 ib0-i0系
        Cni0 = p2cne([pos(1); kts*glv.wie; 0]); qni0 = m2qua(Cni0); %i0-n系
        qnb = qmul(qmul(qni0,qi0ib0),qib0b);
        attk(ki,:) = q2att(qnb)';
        ki = timebar;
    end
    att0 = attk(end,:)';
    resdisp('Initial align attitudes (arcdeg)', att0/glv.deg);
    insplot([attk, (1:length(attk))'*nts]);