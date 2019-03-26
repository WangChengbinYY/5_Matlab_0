function [satPos, satClkCorr, TGD, orbitp] = satPosBatch(transmitTime, eph)
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 30/08/2013, 30/06/2015
global ggps
    s = ggps.ephs;
    
    dt = transmitTime - eph(:,s.Toc);
    idx=dt>302400; dt(idx)=dt(idx)-604800; idx=dt<-302400; dt(idx)=dt(idx)+604800; 
    satClkCorr = eph(:,s.af0) + eph(:,s.af1).*dt+eph(:,s.af2).*dt.*dt;   % �Ӳ������0.1ms������
    
    tk = transmitTime - satClkCorr - eph(:,s.Toe);                    % �黯ʱ��
    idx=dt>302400; tk(idx)=tk(idx)-604800; idx=tk<-302400; tk(idx)=tk(idx)+604800; 

    a   = eph(:,s.sqrtA).^2;      % ������
    n0  = sqrt(ggps.GM./a.^3);    % ƽ��������
    n   = n0 + eph(:,s.Deltan);   % �㶯����

    M   = eph(:,s.M0) + n.*tk;     % ƽ�����
    E   = M;  dE  = 1.0;          % ƫ�����
    while max(abs(dE)) > 1.e-12
        E_old   = E;
        sinE    = sin(E);
        E       = M + eph(:,s.e).*sinE;
        dE      = E-E_old;
    end
    if max(abs(E-eph(:,s.e).*sin(E)-M))>1e-10, error('E error'); end
    cosE = cos(E);
    nu   = atan2(sqrt(1-eph(:,s.e).^2).*sinE, cosE-eph(:,s.e));  % ������
    
    phi  = nu + eph(:,s.omega);   % �����Ǿ�
    sin2phi = sin(2*phi); cos2phi = cos(2*phi);

    u = phi +                             eph(:,s.Cuc).*cos2phi + eph(:,s.Cus).*sin2phi;  % �����Ǿ����
    r = a.*(1-eph(:,s.e).*cosE) +         eph(:,s.Crc).*cos2phi + eph(:,s.Crs).*sin2phi;  % ʧ��������
    i = eph(:,s.i0) + eph(:,s.iDot).*tk + eph(:,s.Cic).*cos2phi + eph(:,s.Cis).*sin2phi;  % �����Ǽ�����

    Omega = eph(:,s.OMEGA0) + (eph(:,s.OMEGADot)-ggps.wie).*tk - ggps.wie*eph(:,s.Toe); % �����㾭��
    
    sinu = sin(u); cosu = cos(u);
    sini = sin(i); cosi = cos(i);
    sinOmega = sin(Omega); cosOmega = cos(Omega);
    x = r.*cosu; y = r.*sinu;  % �����ڹ��ƽ�������꣨ע�⣺X��ָ�������㣬�����ǽ��ص㣩
    satPos = [ x.*cosOmega-y.*cosi.*sinOmega, x.*sinOmega+y.*cosi.*cosOmega, y.*sini ];  % ECEF����
    
    satClkCorr(:,2) = satClkCorr + ggps.F*eph(:,s.e).*eph(:,s.sqrtA).*sinE;   % �Ӳ�����۸�����10ns������������λ��Ӱ��С��
    TGD = eph(:,s.TGD);
    
    orbitp = [dt, tk, i, Omega, nu, u, r];
            