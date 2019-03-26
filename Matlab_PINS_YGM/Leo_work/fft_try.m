xn = IMUdata(:,1)';
% xn = xn - mean(xn);
N = length(xn);
fs = 100;

Xk = fft(xn,N);

Xkabs = abs(Xk);
N2 = fix(N/2);
figure; plot([0:N2]*fs/N,[Xkabs(1),2*Xkabs(2:N2+1)]/N,'-*');