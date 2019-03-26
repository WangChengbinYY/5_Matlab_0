function ahrs = MahonyUpdate(ahrs, gyro, acc, mag, ts)
% AHRS using Mahony method.
%
% Prototype: ahrs = MahonyUpdate(ahrs, wm, vm, mag, ts)
% Inputs: ahrs - AHRS structure array
%        gyro - gyro sample in deg/s
%        acc - accelerometer sample in g (or any unit)
%        mag - magetic output in mGauss (or any unit)
%        ts - sample interval
% Output: ahrs - output AHRS structure array
%
% See also  MahonyInit, QEAHRSUpdate.

% Copyright(c) 2009-2017, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 07/06/2017
    nm = norm(acc);
    if nm>0, acc = acc/nm;  else  acc = [0;0;0];  end
    nm = norm(mag);
    if nm>0, mag = mag/nm;  else  mag = [0;0;0];  end
    bxyz = ahrs.Cnb*mag;
    bxyz(1:2) = [0;norm(bxyz(1:2))];
    wxyz = ahrs.Cnb'*bxyz;
    exyz = cros(ahrs.Cnb(3,:)',acc) + cros(wxyz,mag);
%     exyz = angle3d(ahrs.Cnb(3,:)',acc) + angle3d(wxyz,mag);
    ahrs.exyzInt = ahrs.exyzInt + exyz*ahrs.Ki*ts;
    ahrs.q = qmul(ahrs.q, rv2q((gyro*(pi/180)-ahrs.Kp*exyz-ahrs.exyzInt)*ts));
    ahrs.Cnb = q2mat(ahrs.q);
    ahrs.tk = ahrs.tk + ts;
