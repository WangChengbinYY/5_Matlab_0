function ins = insinit_wcb(att,vn,pos,ts,nn)
% 2017.11.27 王成宾修改
% SINS structure array initialization.
%
% Prototype: ins = insinit(avp0, ts, var1, var2)
% Initialization usages(maybe one of the following methods):
%       ins = insinit(avp0, ts);
%       ins = insinit(avp0, ts, avperr);
%       ins = insinit(qnb0, vn0, pos0, ts);
% Inputs: avp0 - initial avp0 = [att0; vn0; pos0]
%         ts - SIMU sampling interval
%         avperr - avp error setting
% Output: ins - SINS structure array
%
% See also  insupdate, avpset, kfinit.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 22/03/2008, 12/01/2013, 18/03/2014

global glv
    ins = [];
	ins.ts = ts;
    [ins.att, ins.vn, ins.pos] = setvals(att,vn,pos); 
    [ins.qnb, ins.att, ins.Cnb] = attsyn(ins.att);
    ins.eth = earth(ins.pos, ins.vn);
    ins.nn = nn;    %计算周期
    ins.nts = nn*ts;    %计算周期时间
%     ins.fb = [0;0;-glv.g0];
%     ins.wib = zeros(3,1);



