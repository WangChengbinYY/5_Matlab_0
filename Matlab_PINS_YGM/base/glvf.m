function glv1 = glvf(Re, f, wie)
% PSINS Toolbox global variable structure initialization.
%
% Prototype: glv = glvf(Re, f, wie)
% Inputs: Re - the Earth's semi-major axis
%         f - flattening
%         wie - the Earth's angular rate
% Output: glv1 - output global variable structure array
%
% See also  psinsinit.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 14/08/2011, 10/09/2013, 09/03/2014
global glv
    if ~exist('Re', 'var'),  Re = [];  end
    if ~exist('f', 'var'),   f = [];  end
    if ~exist('wie', 'var'), wie = [];  end
    if isempty(Re),  Re = 6378137;  end
    if isempty(f),   f = 1/298.257;  end
    if isempty(wie), wie = 7.2921151467e-5;  end
    glv.Re = Re;                    % the Earth's semi-major axis
    glv.f = f;                      % flattening
    glv.Rp = (1-glv.f)*glv.Re;      % semi-minor axis
    glv.e = sqrt(2*glv.f-glv.f^2); glv.e2 = glv.e^2; % 1st eccentricity
    glv.ep = sqrt(glv.Re^2-glv.Rp^2)/glv.Rp; glv.ep2 = glv.ep^2; % 2nd eccentricity
    glv.wie = wie;                  % the Earth's angular rate 单位弧度
    glv.meru = glv.wie/1000;        % milli earth rate unit 
    glv.g0 = 9.7803267714;          % gravitational force
    glv.mg = 1.0e-3*glv.g0;         % milli g
    glv.ug = 1.0e-6*glv.g0;         % micro g
    glv.mGal = 1.0e-3*0.01;         % milli Gal = 1cm/s^2 约等于 1.0E-6*g0
            %GAL常用于地震工程学中，用来描述地震加速度。gal，称为“伽”或者“盖”，
            %是为纪念第一个重力测量者意大利科学家伽利略（1564~1642）而命名的。
            %重力场的量纲是厘米每二次方秒，规定1cm/s2为1gal；1gal的1/1000为1毫伽，又称为米盖；
            %毫伽的1/1000为1微伽，又称u盖。毫伽是重力勘探的常用单位，由于重力测量精度的提高，用微伽者也较常见。
%待定-------------
    glv.ugpg2 = glv.ug/glv.g0^2;    % ug/g^2
    glv.ws = 1/sqrt(glv.Re/glv.g0); % Schuler frequency
%------------------

    glv.ppm = 1.0e-6;               % parts per million
    glv.deg = pi/180;               % arcdeg
    glv.min = glv.deg/60;           % arcmin
    glv.sec = glv.min/60;           % arcsec
    glv.hur = 3600;                 % time hour (1hur=3600second)

%待定------------- 
    glv.dps = pi/180/1;             % arcdeg / second
    glv.dph = glv.deg/glv.hur;      % arcdeg / hour
    glv.dpss = glv.deg/sqrt(1);     % arcdeg / sqrt(second)
    glv.dpsh = glv.deg/sqrt(glv.hur);  % arcdeg / sqrt(hour)
    glv.dphpsh = glv.dph/sqrt(glv.hur); % (arcdeg/hour) / sqrt(hour)       
    glv.Hz = 1/1;                   % Hertz
    glv.dphpsHz = glv.dph/glv.Hz;   % (arcdeg/hour) / sqrt(Hz)
    glv.ugpsHz = glv.ug/sqrt(glv.Hz);  % ug / sqrt(Hz)
    glv.ugpsh = glv.ug/sqrt(glv.hur); % ug / sqrt(hour)
%        [54,   92,   214,  0,    0    ]/105
    glv.mpsh = 1/sqrt(glv.hur);     % m / sqrt(hour)
    glv.mpspsh = 1/1/sqrt(glv.hur); % (m/s) / sqrt(hour), 1*mpspsh=1667*ugpsHz
    glv.ppmpsh = glv.ppm/sqrt(glv.hur); % ppm / sqrt(hour)
%------------------

    glv.mil = 2*pi/6000;            % mil   密位
    glv.nm = 1853;                  % nautical mile  一海里
    glv.kn = glv.nm/glv.hur;        % knot
           % 节（英文：knot），单位符号kn，是一个专用于航海的速率单位，后延伸至航空方面，相当于船只或飞机每小时所航行的海里数。
           %1节的速度为：每小时1海里（n mile/h）每小时1.852千米（km/h) 每秒0.5144444米（m/s）[1]   
    
    %%
%待定-------------        
    glv.wm_1 = [0,0,0]; glv.vm_1 = [0,0,0];   % the init of previous gyro & acc sample
%------------------
    glv.cs = [                      % coning & sculling compensation coefficients 圆锥划船补偿洗漱
        [2,    0,    0,    0,    0    ]/3
        [9,    27,   0,    0,    0    ]/20
        [250,  525,  650,  1375, 0    ]/504 
        [2315, 4558, 7296, 7834, 15797]/4620 ];
    
 %待定-------------      
    glv.csmax = size(glv.cs,1)+1;  % max subsample number
 %------------------   
    glv.v0 = [0;0;0];    % 3x1 zero-vector
    glv.qI = [1;0;0;0];  % identity quaternion
    glv.I33 = eye(3); glv.o33 = zeros(3);  % identity & zero 3x3 matrices
    glv.pos0 = [34.246048*glv.deg; 108.909664*glv.deg; 380]; % position of INS Lab@NWPU
    glv.eth = []; glv.eth = earth(glv.pos0);
    %%
    [glv.rootpath, glv.datapath, glv.mytestflag] = psinsenvi;
    glv1 = glv;

