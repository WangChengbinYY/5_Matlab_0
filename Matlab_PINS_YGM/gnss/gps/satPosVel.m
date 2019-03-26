function [satPoss, satClkCorrs, satVels] = satPosVel(transmitTime, eph)
% Calculate satellite position(s), clock error(s) and velocity(s) from ephemeris data.
%
% Prototype: [satPoss, satClkCorrs, satVels] = satPosVel(transmitTime, eph) 
% Inputs: transmitTime - satellite signal transmission time
%         eph - ephemeris data
% Outputs: satPoss - satellite positions in ECEF at transmit time
%          satClkCorrs - satellite clock corrections
%
% See also  findEph, lsPos.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 30/08/2013, 30/06/2015
    numOfSatellites = size(transmitTime, 1);
    satPoss         = zeros(3, numOfSatellites); satVels = satPoss;
    satClkCorrs     = zeros(numOfSatellites, 1);
    for k=1:numOfSatellites
        if nargout==3
            [satPoss(:,k), satClkCorrs(k), satVels(:,k)] = satpv(transmitTime(k), eph(k));
        else
            [satPoss(:,k), satClkCorrs(k)] = satpv(transmitTime(k), eph(k));
        end
    end
    
function [satPos, satClkCorr, satVel] = satpv(transmitTime, eph)
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 30/08/2013, 30/06/2015
global ggps
    dt = transmitTime - eph.Toc;
    if dt>302400, dt=dt-604800; elseif dt<-302400, dt=dt+604800;  end
    satClkCorr = eph.af0 + eph.af1*dt+eph.af2*dt*dt - 0*eph.TGD;   % �Ӳ������0.1ms������
    
    tk = transmitTime - satClkCorr - eph.Toe;                    % �黯ʱ��
    if tk>302400, tk=tk-604800; elseif tk<-302400, tk=tk+604800;  end

    a   = eph.sqrtA * eph.sqrtA;   % ������
    n0  = sqrt(ggps.GM / a^3);   % ƽ��������
    n   = n0 + eph.Deltan;   % �㶯����

    M   = eph.M0 + n*tk;   % ƽ�����
    E   = M;  dE  = 1.0;    % ƫ�����
    while abs(dE) > 1.e-12
        E_old   = E;
        sinE    = sin(E);
        E       = M + eph.e*sinE;
        dE      = E-E_old;
    end
    if abs(E-eph.e*sin(E)-M)>1e-10, error('E error'); end
    cosE = cos(E);
    nu   = atan2(sqrt(1-eph.e^2)*sinE, cosE-eph.e);  % ������
    
    phi  = nu + eph.omega;   % �����Ǿ�
    sin2phi = sin(2*phi); cos2phi = cos(2*phi);

    u = phi +                  eph.Cuc*cos2phi + eph.Cus*sin2phi;  % �����Ǿ����
    r = a*(1-eph.e*cosE) +     eph.Crc*cos2phi + eph.Crs*sin2phi;  % ʧ��������
    i = eph.i0 + eph.iDot*tk + eph.Cic*cos2phi + eph.Cis*sin2phi;  % �����Ǽ�����

    Omega = eph.OMEGA0 + (eph.OMEGADot-ggps.wie)*tk - ggps.wie*eph.Toe; % �����㾭��
    
    sinu = sin(u); cosu = cos(u);
    sini = sin(i); cosi = cos(i);
    sinOmega = sin(Omega); cosOmega = cos(Omega);
    x = r*cosu; y = r*sinu;  % �����ڹ��ƽ�������꣨ע�⣺X��ָ�������㣬�����ǽ��ص㣩
    satPos = [ x*cosOmega - y*cosi*sinOmega; x*sinOmega + y*cosi*cosOmega; y*sini ];  % ECEF����
    
    satClkCorr = satClkCorr + ggps.F*eph.e*eph.sqrtA*sinE;   % �Ӳ�����۸�����10ns������������λ��Ӱ��С��
            
    if nargout==3      % by Ling Yang 2012/5/24 at UNSW
        E_dot = n/(1-eph.e*cosE);
        v_dot = E_dot*(sqrt(1-eph.e^2))/(1-eph.e*cosE);
        r_dot = a*eph.e*sinE*E_dot + (2*eph.Crs*cos2phi-2*eph.Crc*sin2phi)*v_dot;
        u_dot = (1+2*eph.Cus*cos2phi-2*eph.Cuc*sin2phi)*v_dot;
        W_dot = eph.OMEGADot - ggps.wie;
        i_dot = 2*(eph.Cis*cos2phi-eph.Cic*sin2phi)*v_dot + eph.iDot;
        y_dot = r_dot*sinu + r*cosu*u_dot;
        x_dot = r_dot*cosu - r*sinu*u_dot;
        satVel = [...
        	x_dot*cosOmega - y_dot*cosi*sinOmega + y*sini*sinOmega*i_dot - satPos(2)*W_dot; ...
        	x_dot*sinOmega + y_dot*cosi*cosOmega - y*sini*cosOmega*i_dot + satPos(1)*W_dot; ...
        	y_dot*sini + y*cosi*i_dot ];
    end