function ahrs = MahonyInit(tau)
% The initialization of AHRS using Mahony method. The transfer function is
% as follow:
%         s*Kp+Ki                s
%  dq = ----------- * fb  +  ----------- * wb
%       s^2+s*Kp+Ki          s^2+s*Kp+Ki
%
% Prototype: ahrs = MahonyInit(tau)
% Inputs: tau - time constant, Delta(s) = s^2+Kp*s+Ki
%                                       = s^2+2*beta*s+beta^2
%                              where beta = 1/tau
% Output: ahrs - output AHRS structure array
%
% See also  MahonyUpdate, QEAHRSInit.

% Copyright(c) 2009-2017, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 07/06/2017
    ahrs.q = [1;0;0;0];  ahrs.Cnb = q2mat(ahrs.q);
    ahrs.exyzInt = [0;0;0];
    if ~exist('tau', 'var'), tau = 4; end
    beta = 2.146/tau;
    ahrs.Kp = 2*beta; ahrs.Ki = beta^2;
    ahrs.tk = 0;

