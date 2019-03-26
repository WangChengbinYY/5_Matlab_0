function [obs, head] = rnx210o(oFileName, nRecord)
% A simple version of read RINEX-format (Ver2.10) observation data file.
%
% Prototype: [obs, head] = rnx210o(oFileName, nRecord)
% Inputs: oFileName - RINEX-format GPS observation file name
%         nRecord - first n records to be readed
% Outputs: obss - observation records in struct
%          obsTpIdx - observation record index = observation time array [tp]
%          head - obs file header information
%
% See also  rnx210n, obsEphLink.

% Originated by Yangling 2008/3/25, Tongji Unv.,China
% Modified by Yangongmin 16/09/2013, 05/07/2015, NWPU
%% read header
	fid = fopen(oFileName, 'rt');
    s = fgetl(fid);
    if ~strfind(s,'OBSERVATION DATA')
        error('can not find the beginning of header!\r\n');
    end
    for lines=2:100
        s = fgetl(fid); s(length(s)+1:81) = ' ';
        if strcmp(s(61:71),'MARKER NAME')  %��վ���   
            head.siteName = s(1:4);    
        elseif strcmp(s(61:79),'APPROX POSITION XYZ') %��վ���������
            head.xyz0 = str2num(s(1:42));
        elseif strcmp(s(61:80),'ANTENNA: DELTA H/E/N') %�������ĸ��� 
            head.antenna = str2num(s(1:42));
        elseif strcmp(s(61:79),'# / TYPES OF OBSERV') %�۲�ֵ����
            head.obsNum = str2num(s(1:6));
            otype.WN=1; otype.tp=2; otype.PRN=3;
            for k=1:head.obsNum
                head.obsName(k).Name = s(k*6+5:k*6+6);
                otype = setfield(otype, s(k*6+5:k*6+6), 3+k);
            end
            head.otype = otype;
        elseif strcmp(s(61:68),'INTERVAL') %��Ԫ���
            head.interval=str2num(s(5:11));
        elseif strcmp(s(61:77),'TIME OF FIRST OBS') %��һ����Ԫʱ�� 
            head.timeOB=cal2gpst(str2num(s(1:42)));
        elseif strcmp(s(61:76),'TIME OF LAST OBS') %���һ����Ԫʱ�� 
            head.timeOE=transTime(s(1:42));
        elseif strcmp(s(61:73),'END OF HEADER') %ͷ�ļ����� 
            break;
        end 
    end
%% read obs
    if nargin<2,  nRecord = 20000; end
    obs = [];  idx = 1;
    for k=1:nRecord
        obsi = obs_read(fid, head);
        if feof(fid), break;  end
        if length(obs)<idx+obsi.satNum,  obs=[obs; zeros(1000,head.obsNum+3)]; end
        for kk=1:obsi.satNum
            obs(idx,:) = [obsi.wn,obsi.tp,obsi.PRNs(kk),obsi.obs(kk,:)];
            idx = idx + 1;
        end
    end
    obs(idx:end,:) = [];
    fclose(fid);
  

function obs = obs_read(fid, head)
	while 1   %ÿ����Ԫ�ĵ�һ�����ݣ�ʱ��͹۲⵽�����Ǻ�
        s = fgetl(fid); s(length(s)+1:81) = ' ';
        if isempty(str2num(s(1:3)))
            continue;
        end
        [tow, wn] = cal2gpst(str2num(s(1:26))); %��Ԫʱ��
        obs.wn = wn;
        obs.tp = tow;
        obs.flag = str2num(s(27:29)); %Epoch flag
        if isempty(obs.flag)
            obs.flag = 7;
        end
        if obs.flag==4 
            s = fgetl(fid); s(length(s)+1:81) = ' '; %һ��Сʱ�����ݽ�������һ��Сʱ�����ݿ�ʼ������Ϊ˵������
            continue;
        else
            break;
        end
	end
    ss = s;
    obs.satNum = str2num(s(30:32));  %����Ԫ�۲⵽��������
    s(33:3:end) = ' ';
    obs.PRNs = str2num(s(33:end));
    if length(obs.PRNs)<obs.satNum
        s1 = fgetl(fid); s1(length(s1)+1:81) = ' ';
        s = [s(1:33+3*length(obs.PRNs)-1), s1(33:end)];
        ss = [ss, s1];
        s(33:3:end) = ' ';
        obs.PRNs = str2num(s(33:end));
    end

    obs.dot = str2num(s(69:80)); %���ջ��Ӳ�d
    if isempty(obs.dot)
        obs.dot=0;
    end

    obs.obs = zeros(obs.satNum,head.obsNum);
    for i=1:obs.satNum  %ÿ����Ԫ�Ĺ۲����ݣ������Ǻ��Ⱥ�˳����д洢
        if feof(fid)
            break;
            obs.flag = 9;
        end
        s = fgetl(fid); s(length(s)+1:81) = ' ';
        if head.obsNum>5  %�ж�һ�����ǵĹ۲������Ƿ�����д洢�����Ϊ���У����ٶ���һ��
            s1 = fgetl(fid); s1(length(s1)+1:81) = ' ';
            s = [s, s1]; 
        end
%         obs.obs(i,:) = [str2num(s),0];
        obs.obs(i,:) = s2n(s, head.obsNum);
    end

    if obs.satNum<4
        obs.flag = 8;    
    end

	obs.satNum = length(find(ss=='G'));
    obs.PRNs = obs.PRNs(1:obs.satNum);
    obs.obs = obs.obs(1:obs.satNum,:);
    
function num = s2n(str, n)
    num = zeros(1, n);
    for k = 1:n
        kk = (k-1)*16+2;
        numi = str2num(str(kk:kk+13));
        if isempty(numi)
            numi = 0;
        end
        num(k) = numi;
    end