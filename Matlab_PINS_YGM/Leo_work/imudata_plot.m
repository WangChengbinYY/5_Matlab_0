function imudata_plot(IMUdata)
% ����imu���ݵ�ͼ��

figure;
subplot(231);
plot(attfirst(:,4),attfirst(:,1)),
title('����'),xlabel('����'),ylabel('����/��');
subplot(232);
plot(attfirst(:,4),attfirst(:,2)),
title('���'),xlabel('����'),ylabel('���/��');
subplot(233);
plot(attfirst(:,4),attfirst(:,3)),
title('����'),xlabel('����'),ylabel('����/��');

subplot(234);
plot(attsecond(:,4),attsecond(:,1)),
title('����'),xlabel('����'),ylabel('����/��');
subplot(235);
plot(attsecond(:,4),attsecond(:,2)),
title('���'),xlabel('����'),ylabel('���/��');
subplot(236);
plot(attsecond(:,4),attsecond(:,3)),
title('����'),xlabel('����'),ylabel('����/��');






%%  ����ԭʼ����
t = [1:1:L]*ts;     %������Ϊʱ��
figure;
subplot(231);
plot(t,IMUdata(:,1)),axis([0 L*ts minmax(IMUdata(:,1)')]),
title('X������'),xlabel('t/s'),ylabel('����/s');
subplot(232);
plot(t,IMUdata(:,2)),axis([0 L*ts minmax(IMUdata(:,2)')]),
title('Y������'),xlabel('t/s'),ylabel('����/s');
subplot(233);
plot(t,IMUdata(:,3)),axis([0 L*ts minmax(IMUdata(:,3)')]),
title('Z������'),xlabel('t/s'),ylabel('����/s');
subplot(234);
plot(t,IMUdata(:,4)),axis([0 L*ts minmax(IMUdata(:,4)')]),
title('X����ٶȼ�'),xlabel('t/s'),ylabel('m/s2');
subplot(235);
plot(t,IMUdata(:,5)),axis([0 L*ts minmax(IMUdata(:,5)')]),
title('Y����ٶȼ�'),xlabel('t/s'),ylabel('m/s2');
subplot(236);
plot(t,IMUdata(:,6)),axis([0 L*ts minmax(IMUdata(:,6)')]),
title('Z����ٶȼ�'),xlabel('t/s'),ylabel('m/s2');

%%  ��ԭʼ����ֱ�ӻ��� ���� �Ӽ�
IMUdata_Integrate = zeros(L,6);
IMUdata_Integrate(1,1:5) = IMUdata(1,1:5)*ts;
IMUdata_Integrate(1,6) = (IMUdata(1,6)-glv.g0)*ts;
for i=2:L
    IMUdata_Integrate(i,1:5) = IMUdata_Integrate(i-1,1:5) +  IMUdata(i,1:5)*ts;    
    IMUdata_Integrate(i,6) = IMUdata_Integrate(i-1,6) +  (IMUdata(i,6)-glv.g0)*ts;    
end;
    IMUdata_Integrate(:,1:3) = IMUdata_Integrate(:,1:3)/glv.deg;
    
figure;
subplot(231);
plot(t,IMUdata_Integrate(:,1)),axis([0 L*ts minmax(IMUdata_Integrate(:,1)')]),
title('X������ʱ�����'),xlabel('t/s'),ylabel('��');
subplot(232);
plot(t,IMUdata_Integrate(:,2)),axis([0 L*ts minmax(IMUdata_Integrate(:,2)')]),
title('Y������ʱ�����'),xlabel('t/s'),ylabel('��');
subplot(233);
plot(t,IMUdata_Integrate(:,3)),axis([0 L*ts minmax(IMUdata_Integrate(:,3)')]),
title('Z������ʱ�����'),xlabel('t/s'),ylabel('��');
subplot(234);
plot(t,IMUdata_Integrate(:,4)),axis([0 L*ts minmax(IMUdata_Integrate(:,4)')]),
title('X����ٶȼ�ʱ�����'),xlabel('t/s'),ylabel('m');
subplot(235);
plot(t,IMUdata_Integrate(:,5)),axis([0 L*ts minmax(IMUdata_Integrate(:,5)')]),
title('Y����ٶȼ�ʱ�����'),xlabel('t/s'),ylabel('m');
subplot(236);
plot(t,IMUdata_Integrate(:,6)),axis([0 L*ts minmax(IMUdata_Integrate(:,6)')]),
title('Z����ٶȼ�ʱ�����'),xlabel('t/s'),ylabel('m');
