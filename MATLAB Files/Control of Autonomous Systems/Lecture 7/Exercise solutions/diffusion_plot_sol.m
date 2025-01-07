% plot results of diffusion
% Jerome Jouffroy, November 2023

x_trace = out.x;
[nt,nxi] = size(x_trace);

t = 0:nt-1;
ix = 1:nxi;

figure, surface(ix,t,x_trace,'FaceColor','r','FaceAlpha',0.7)
xlabel('State components','FontSize',24)
ylabel('Time','FontSize',24)
view(-225,30)
grid on

hold on
x_hat_trace = out.x_hat;
surface(ix,t,x_hat_trace,'FaceColor','g')
hold off
legend('state','observer','FontSize',24,'Location','northeast')