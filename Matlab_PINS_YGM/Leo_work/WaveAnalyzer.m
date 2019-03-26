%waveletAnalyzer

gyro_x = IMUdata(:,4);
gyro_y = IMUdata(:,4);
gyro_z = IMUdata(:,4);

IMUdataNew = IMUdata;
IMUdataNew(:,4) = gyrox_new';
IMUdataNew(:,5) = gyroy_new';
IMUdataNew(:,6) = gyroz_new';

L = length(imudata);
G = zeros(L,1);
for i=1:L
    G(i,1) = sqrt(imudata(i,4)^2+imudata(i,5)^2+imudata(i,6)^2);
end
    