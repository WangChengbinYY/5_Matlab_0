function magplot(mag)
% 3-axis geomagnetic plot.
%
% Prototype: magplot(mag)
% Input: mag - 3-axis geomagnetic in milli-Gauss
%          
% See also  imuplot, insplot, inserrplot, kfplot, gpsplot.

% Copyright(c) 2009-2017, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/03/2017
    myfigure;
    plot(mag(:,end), [mag(:,1:3),normv(mag(:,1:3))]), xygo('\itt/ \ims', 'Mag / mGauss');
