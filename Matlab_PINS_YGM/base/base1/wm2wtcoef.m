function coef = wm2wtcoef(ts, n, N)
% Determine the transform matrix form angular increment to angular rate.
%
% Prototype: coef = wm2wtcoef(ts, n, N)
% Inputs: ts - sample interval
%         n - current sub-sample number
%         N - total sub-sample number, then p=N-n is previous sub-sample
%             nmuber
% Output: coef - angular rate coefficient matrix, such that wt = wm'*coef
%
% See also  qpicard, dcmtaylor.

% Copyright(c) 2009-2017, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 04/02/2017
    if nargin<3, N = n; end
    G = zeros(N); p = N - n;
    for k=1:N
        k1 = N-k+1;
        for j=(-p+1):n
            G(k,j+p) = ((j*ts)^k1-((j-1)*ts)^k1) / k1;
        end
    end
    coef = G^-1;