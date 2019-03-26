function dop = DOP(A)
% Calculate GPS positioning DOP values.
%
% See also  lsPos.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 30/08/2013, 30/06/2015
    Q = inv(A'*A);
    dop = sqrt([...
        Q(1,1) + Q(2,2) + Q(3,3) + Q(4,4);       % GDOP
        Q(1,1) + Q(2,2) + Q(3,3);                % PDOP
        Q(1,1) + Q(2,2);                         % HDOP
        Q(3,3);                                  % VDOP
        Q(4,4);                                  % TDOP
        ]);
        
