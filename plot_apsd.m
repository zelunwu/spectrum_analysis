function plot_apsd(freq,apsd,alpha,m_segment,x_frac,y_frac)


% Calculate the confidence level.
df = 2*m_segment;

x1 = chi2inv(1-alpha/2,df);
x2 = chi2inv(alpha/2,df);

sign_min = log10(df/x1);
sign_max = log10(df/x2);


%  plot figure;
% figure;
h = loglog(freq,apsd,'linewidth',2);
hold on

xlog_range = log10(get(gca,'xlim'));
ylog_range = log10(get(gca,'ylim'));

sign_point_loc_x = min(xlog_range) + range(xlog_range)*x_frac;
sign_point_loc_y = min(ylog_range) + range(ylog_range)*y_frac;

sign_point_y_min = sign_point_loc_y + sign_min;
sign_point_y_max = sign_point_loc_y + sign_max;

plot(10.^(ones(3,1)*sign_point_loc_x),10.^([sign_point_y_min,sign_point_loc_y,sign_point_y_max]),'o-','linewidth',2,'Color',h.Color);
text(10^(sign_point_loc_x+range(xlog_range)*0.02),10^sign_point_loc_y,[num2str((1-alpha)*100),'%, df = ',num2str(df)],'fontsize',15);
% hold off
% xlabel('Frequency (Hz)');
% ylabel('Power Spectrum Density');

set(gca,'fontsize',15);
% hold off

end