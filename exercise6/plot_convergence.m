error = load('error.txt');
dt = load('dt.txt');

N = length(error);
rate = polyfit(log(dt(2:N)), log(error(2:N)), 1);
rate = rate(1);
disp(rate);

loglog(dt, error, '-o')
hold on
loglog(dt, error(N-1)*dt.^rate/(dt(N-1).^rate), '--')
xlabel('\Delta t')
ylabel('Error')
axis([dt(N) dt(1) -inf inf]);
label = strcat('\Delta t^{', num2str(rate,3), '}');
legend('|u(T)-U_N|', label, 'Location', 'northwest')
grid on
