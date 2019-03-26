function [satPos, satClkCorr, TGD, orbitp] = bdsatPosBatch(transmitTime, eph)
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/07/2015
global gbd
    s = gbd.ephs;
    
    dt = transmitTime - eph(:,s.Toc);
    idx=dt>302400; dt(idx)=dt(idx)-604800; idx=dt<-302400; dt(idx)=dt(idx)+604800; 
    satClkCorr = eph(:,s.af0) + eph(:,s.af1).*dt+eph(:,s.af2).*dt.*dt;   % �Ӳ������0.1ms������
    
    tk = transmitTime - satClkCorr - eph(:,s.Toe);                    % �黯ʱ��
    idx=dt>302400; tk(idx)=tk(idx)-604800; idx=tk<-302400; tk(idx)=tk(idx)+604800; 

    a   = eph(:,s.sqrtA).^2;      % ������
    n0  = sqrt(gbd.GM./a.^3);    % ƽ��������
    n   = n0 + eph(:,s.Deltan);   % �㶯����

    M   = eph(:,s.M0) + n.*tk;     % ƽ�����
    E   = M;  dE  = 1.0;          % ƫ�����
    while max(abs(dE)) > 1.e-12
        E_old   = E;
        sinE    = sin(E);
        E       = M + eph(:,s.e).*sinE;
        dE      = E-E_old;
    end
    cosE = cos(E);
    nu   = atan2(sqrt(1-eph(:,s.e).^2).*sinE, cosE-eph(:,s.e));  % ������
    
    phi  = nu + eph(:,s.omega);   % �����Ǿ�
    sin2phi = sin(2*phi); cos2phi = cos(2*phi);

    u = phi +                             eph(:,s.Cuc).*cos2phi + eph(:,s.Cus).*sin2phi;  % �����Ǿ����
    r = a.*(1-eph(:,s.e).*cosE) +         eph(:,s.Crc).*cos2phi + eph(:,s.Crs).*sin2phi;  % ʧ��������
    i = eph(:,s.i0) + eph(:,s.iDot).*tk + eph(:,s.Cic).*cos2phi + eph(:,s.Cis).*sin2phi;  % �����Ǽ�����

    sinu = sin(u); cosu = cos(u);
    sini = sin(i); cosi = cos(i);
    x = r.*cosu; y = r.*sinu;  % �����ڹ��ƽ�������꣨ע�⣺X��ָ�������㣬�����ǽ��ص㣩
    % satellite is GEO
    isGEO = eph(:,s.PRN)<=505 & eph(:,s.PRN)>=500;
    if ~isempty(tk(isGEO))
        Omega = eph(isGEO,s.OMEGA0) + eph(isGEO,s.OMEGADot).*tk(isGEO) - gbd.wie*eph(isGEO,s.Toe);
        sinOmega = sin(Omega); cosOmega = cos(Omega);
        xyz = [ x(isGEO).*cosOmega-y(isGEO).*cosi(isGEO).*sinOmega, x(isGEO).*sinOmega+y(isGEO).*cosi(isGEO).*cosOmega, y(isGEO).*sini(isGEO) ];
        wtk = gbd.wie*tk(isGEO);
        swtk = sin(wtk);      cwtk = cos(wtk);
        s5 = sin(-5*pi/180);  c5 = cos(-5*pi/180);
        satPos = zeros(length(tk),3);
        satPos(isGEO,1) =  cwtk.*xyz(:,1) + swtk.*c5.*xyz(:,2) + swtk.*s5.*xyz(:,3);
        satPos(isGEO,2) = -swtk.*xyz(:,1) + cwtk.*c5.*xyz(:,2) + cwtk.*s5.*xyz(:,3);
        satPos(isGEO,3) =               0 -       s5.*xyz(:,2) +       c5.*xyz(:,3);
    end
    % satellite is IGSO/MEO
    noGEO = ~isGEO;
    if ~isempty(tk(noGEO))
        Omega = eph(noGEO,s.OMEGA0) + (eph(noGEO,s.OMEGADot)-gbd.wie).*tk(noGEO) - gbd.wie*eph(noGEO,s.Toe); % �����㾭��
        sinOmega = sin(Omega); cosOmega = cos(Omega);
        satPos(noGEO,:) = [ x(noGEO).*cosOmega-y(noGEO).*cosi(noGEO).*sinOmega, x(noGEO).*sinOmega+y(noGEO).*cosi(noGEO).*cosOmega, y(noGEO).*sini(noGEO) ];  % ECEF����
    end
    
    satClkCorr(:,2) = satClkCorr + gbd.F*eph(:,s.e).*eph(:,s.sqrtA).*sinE;   % �Ӳ�����۸�����10ns������������λ��Ӱ��С��
    TGD = [eph(:,s.TGD1), eph(:,s.TGD2), (eph(:,s.TGD2)-gbd.kf*eph(:,s.TGD1))/(1-gbd.kf)];

	orbitp = [dt, tk, i, nu, u, r];
