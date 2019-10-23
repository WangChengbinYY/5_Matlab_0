function imu = imuadderr(imu, imuerr, ts)
% SIMU outputs adding errors simulation, denoted as:
%    imu = K*imu + drift error.
%
% Prototype: imu = imuadderr(imu, imuerr, ts)
% Inputs: imu - raw SIMU data
%         imuerr - SIMU error struture array
% Output: imu - output SIMU data added errors
%
% See also  imuerrset, avpseterr, trjsimu, insupdate.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 11/09/2013, 06/03/2014, 21/07/2015
    if nargin<3,  ts = imu(2,7)-imu(1,7);  end  % the last column implies sampling interval
    [m, n] = size(imu); sts = sqrt(ts);
    drift = [ ts*imuerr.eb(1) + sts*imuerr.web(1)*randn(m,1), ...
              ts*imuerr.eb(2) + sts*imuerr.web(2)*randn(m,1), ...
              ts*imuerr.eb(3) + sts*imuerr.web(3)*randn(m,1), ...
              ts*imuerr.db(1) + sts*imuerr.wdb(1)*randn(m,1), ...
              ts*imuerr.db(2) + sts*imuerr.wdb(2)*randn(m,1), ...
              ts*imuerr.db(3) + sts*imuerr.wdb(3)*randn(m,1) ];
%     if imuerr.sqg(1)>0
%         mvg = markov1(imuerr.sqg, imuerr.taug, ts, m);
%         drift(:,1:3) = drift(:,1:3) + mvg*ts;
%     end
%     if imuerr.sqa(1)>0
%         mva = markov1(imuerr.sqa, imuerr.taua, ts, m);
%         drift(:,4:6) = drift(:,4:6) + mva*ts;
%     end
     if min(abs(imuerr.sqg))>0
        mvg = markov1(imuerr.sqg*sqrt(imuerr.taug/2), imuerr.taug, ts, m);   % q = 2*sigma.^2.*beta
        drift(:,1:3) = drift(:,1:3) + mvg*ts;
     end
     if min(abs(imuerr.sqa))>0
        mva = markov1(imuerr.sqa*sqrt(imuerr.taua/2), imuerr.taua, ts, m);
        drift(:,4:6) = drift(:,4:6) + mva*ts;
     end

    if isfield(imuerr, 'dKg')
        Kg = eye(3)+imuerr.dKg; Ka = eye(3)+imuerr.dKa;
        imu(:,1:6) = [imu(:,1:3)*Kg', imu(:,4:6)*Ka'];
    end
    imu(:,1:6) = imu(:,1:6) + drift;
