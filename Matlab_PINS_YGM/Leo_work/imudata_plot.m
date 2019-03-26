function imudata_plot(IMUdata)
% 绘制imu数据的图形

figure;
subplot(231);
plot(attfirst(:,4),attfirst(:,1)),
title('俯仰'),xlabel('分钟'),ylabel('俯仰/度');
subplot(232);
plot(attfirst(:,4),attfirst(:,2)),
title('横滚'),xlabel('分钟'),ylabel('横滚/度');
subplot(233);
plot(attfirst(:,4),attfirst(:,3)),
title('航向'),xlabel('分钟'),ylabel('航向/度');

subplot(234);
plot(attsecond(:,4),attsecond(:,1)),
title('俯仰'),xlabel('分钟'),ylabel('俯仰/度');
subplot(235);
plot(attsecond(:,4),attsecond(:,2)),
title('横滚'),xlabel('分钟'),ylabel('横滚/度');
subplot(236);
plot(attsecond(:,4),attsecond(:,3)),
title('航向'),xlabel('分钟'),ylabel('航向/度');






%%  绘制原始数据
t = [1:1:L]*ts;     %横坐标为时间
figure;
subplot(231);
plot(t,IMUdata(:,1)),axis([0 L*ts minmax(IMUdata(:,1)')]),
title('X轴陀螺'),xlabel('t/s'),ylabel('弧度/s');
subplot(232);
plot(t,IMUdata(:,2)),axis([0 L*ts minmax(IMUdata(:,2)')]),
title('Y轴陀螺'),xlabel('t/s'),ylabel('弧度/s');
subplot(233);
plot(t,IMUdata(:,3)),axis([0 L*ts minmax(IMUdata(:,3)')]),
title('Z轴陀螺'),xlabel('t/s'),ylabel('弧度/s');
subplot(234);
plot(t,IMUdata(:,4)),axis([0 L*ts minmax(IMUdata(:,4)')]),
title('X轴加速度计'),xlabel('t/s'),ylabel('m/s2');
subplot(235);
plot(t,IMUdata(:,5)),axis([0 L*ts minmax(IMUdata(:,5)')]),
title('Y轴加速度计'),xlabel('t/s'),ylabel('m/s2');
subplot(236);
plot(t,IMUdata(:,6)),axis([0 L*ts minmax(IMUdata(:,6)')]),
title('Z轴加速度计'),xlabel('t/s'),ylabel('m/s2');

%%  对原始数据直接积分 陀螺 加计
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
title('X轴陀螺时间积分'),xlabel('t/s'),ylabel('度');
subplot(232);
plot(t,IMUdata_Integrate(:,2)),axis([0 L*ts minmax(IMUdata_Integrate(:,2)')]),
title('Y轴陀螺时间积分'),xlabel('t/s'),ylabel('度');
subplot(233);
plot(t,IMUdata_Integrate(:,3)),axis([0 L*ts minmax(IMUdata_Integrate(:,3)')]),
title('Z轴陀螺时间积分'),xlabel('t/s'),ylabel('度');
subplot(234);
plot(t,IMUdata_Integrate(:,4)),axis([0 L*ts minmax(IMUdata_Integrate(:,4)')]),
title('X轴加速度计时间积分'),xlabel('t/s'),ylabel('m');
subplot(235);
plot(t,IMUdata_Integrate(:,5)),axis([0 L*ts minmax(IMUdata_Integrate(:,5)')]),
title('Y轴加速度计时间积分'),xlabel('t/s'),ylabel('m');
subplot(236);
plot(t,IMUdata_Integrate(:,6)),axis([0 L*ts minmax(IMUdata_Integrate(:,6)')]),
title('Z轴加速度计时间积分'),xlabel('t/s'),ylabel('m');
