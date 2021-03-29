load et_noWalkers
et = et_noWalkers(:,1);

max_lag = 100;
N = length(et);
avg_et = mean(et);
r = zeros(20,1);
for k=1:max_lag
    top = 0;
    bottom = 0;
    for i=1:(N-k)
        top = top + (et(i)-avg_et)*(et(i+k)-avg_et);
    end
    for i=1:N
        bottom = bottom + (et(i) - avg_et)^2;
    end
    r(k) = top/bottom;
end

bar(1:max_lag,r)
xlabel('Lag','Interpreter','latex')
ylabel('Energy autocorrelation','Interpreter','latex')
%title('Autocorrelation of the energy')

% [ac,lag] = xcorr(et,length(et),'coeff');
% plot(lag/length(et),ac)
% plot(1:length(ac),ac)