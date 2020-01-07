function plot_cpsd(freq,coherence,alpha,m_segment)


df = 2*m_segment;
x = chi2inv(1-alpha,2);

coherence_alpha = x/df;

% figure;
plot(freq,coherence,'linewidth',2);
hold on
plot(freq,ones(size(freq))*coherence_alpha,'linewidth',2);
hold off
ylim([0,1]);
legend({'Coherence',['95% Sinificance Level, df = ',num2str(df)]},'fontsize',15);
ylabel('Coherence');
title('Magnitude Square Coherence');
set(gca,'fontsize',15);
end