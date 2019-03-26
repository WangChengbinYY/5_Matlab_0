function imuparameter = Leo_imuparameter(eb, db, web, wdb, sqrtR0G, TauG, sqrtR0A, TauA, dKGii, dKAii, dKGij, dKAij, KA2)
% SIMU errors setting, including gyro & acc bias, noise and installation errors, etc.
%
% Prototype: imuparameter = Leo_imuparameter(eb, db, web, wdb, sqrtR0G, TauG, sqrtR0A, TauA, dKGii, dKAii, dKGij, dKAij)
% Inputs: including infomation as follows
%     eb - gyro constant bias (deg/h)
%     db - acc constant bias (ug)
%     web - angular random walk (deg/sqrt(h))
%     wdb - velocity random walk (ug/sqrt(Hz))
%     sqrtR0G,TauG - gyro correlated bias, sqrtR0G in deg/h and TauG in s
%     sqrtR0A,TauA - acc correlated bias, sqrtR0A in ug and TauA in s
%     dKGii - gyro scale factor error (ppm)
%     dKAii - acc scale factor error (ppm)
%     dKGij - gyro installation error (arcsec)
%     dKAij - acc installation error (arcsec)
%     KA2 - acc quadratic coefficient (ug/g^2)
%     where, 
%             |dKGii(1) dKGij(4) dKGij(5)|         |dKAii(1) 0        0       |
%       dKg = |dKGij(1) dKGii(2) dKGij(6)| , dKa = |dKAij(1) dKAii(2) 0       |
%             |dKGij(2) dKGij(3) dKGii(3)|         |dKAij(2) dKAij(3) dKAii(3)|
% Output: imuparameter - SIMU error structure array
%
% Example:
%     For inertial grade SIMU, typical errors are:
%       eb=0.01dph, db=50ug, web=0.001dpsh, wdb=10ugpsHz
%       scale factor error=10ppm, askew installation error=10arcsec
%       sqrtR0G=0.001dph, taug=1000s, sqrtR0A=10ug, taug=1000s
%    then call this funcion by
%       imuparameter = Leo_imuparameter(0.01,100,0.001,10, 0.001,1000,10,1000, 10,10,10,10);
%
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China

% 10/20/2017 Modified by Chengbin Wang, Tsinghua University
%
%

%% 
%--------------------------------------------------------
%判断IMU系数误差至少输入 陀螺、加计的零偏和零偏稳定性。
%---------------------------------------------------------
global glv
if nargin <4
   msgbox('请输入陀螺、加计的零偏和零偏稳定性!');  
   return;
end
%建立IMU系统参数结构体，num代表参数个数，一般为4
o31 = zeros(3,1); o33 = zeros(3);
imuparameter = struct('num',0,'eb',o31, 'db',o31, 'web',o31, 'wdb',o31,...
    'sqg',o31, 'taug',inf(3,1), 'sqa',o31, 'taua',inf(3,1), 'dKg',o33, 'dKa',o33, 'dKga',zeros(15,1),'KA2',o31); 
switch nargin
    case 4
        imuparameter.num = 4;
        imuparameter.eb(1:3) = eb*glv.dph;   
        imuparameter.web(1:3) = web*glv.dpsh;
        imuparameter.db(1:3) = db*glv.ug;    
        imuparameter.wdb(1:3) = wdb*glv.ugpsHz;
    otherwise
       %% correlated bias
        if exist('sqrtR0G', 'var')
            if TauG(1)==inf, imuparameter.sqg(1:3) = sqrtR0G*glv.dphpsh;   % algular rate random walk !!!
            elseif TauG(1)>0, imuparameter.sqg(1:3) = sqrtR0G*glv.dph.*sqrt(2./TauG); imuparameter.taug(1:3) = TauG; % Markov process
            end
            imuparameter.num = imuparameter.num+2;
        end
        if exist('sqrtR0A', 'var')
            if TauA(1)==inf, imuparameter.sqa(1:3) = sqrtR0A*glv.ugpsh;   % specific force random walk !!!
            elseif TauA(1)>0, imuparameter.sqa(1:3) = sqrtR0A*glv.ug.*sqrt(2./TauA); imuparameter.taua(1:3) = TauA; % Markov process
            end
            imuparameter.num = imuparameter.num+2;
        end
        %% scale factor error
        if exist('dKGii', 'var')
            imuparameter.dKg = setdiag(imuparameter.dKg, dKGii*glv.ppm);
            imuparameter.num = imuparameter.num+1;
        end
        if exist('dKAii', 'var')
            imuparameter.dKa = setdiag(imuparameter.dKa, dKAii*glv.ppm);
            imuparameter.num = imuparameter.num+1;
        end
        %% installation angle error
        if exist('dKGij', 'var')
            dKGij = ones(6,1).*dKGij*glv.sec;
            imuparameter.dKg(2,1) = dKGij(1); imuparameter.dKg(3,1) = dKGij(2); imuparameter.dKg(3,2) = dKGij(3); 
            imuparameter.dKg(1,2) = dKGij(4); imuparameter.dKg(1,3) = dKGij(5); imuparameter.dKg(2,3) = dKGij(6);
            imuparameter.num = imuparameter.num+1;
        end
        if exist('dKAij', 'var')
            dKAij = ones(3,1).*dKAij*glv.sec;
            imuparameter.dKa(2,1) = dKAij(1); imuparameter.dKa(3,1) = dKAij(2); imuparameter.dKa(3,2) = dKAij(3); 
            imuparameter.num = imuparameter.num+1;
        end
        imuparameter.dKga = [imuparameter.dKg(:,1); imuparameter.dKg(:,2);   imuparameter.dKg(:,3);
                       imuparameter.dKa(:,1); imuparameter.dKa(2:3,2); imuparameter.dKa(3,3)];
        if exist('KA2', 'var')
            imuparameter.KA2(1:3) = KA2*glv.ugpg2; 
            imuparameter.num = imuparameter.num+1;
        end         
end

