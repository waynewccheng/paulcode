% Q: Ocean Optics QE65-PRO measurement data validation

figure('Unit','normalized','Position',[0 0 1 1])

k = 1;
for wl = 380:10:780
    
    subplot(6,7,k)
    fn = sprintf('%dnm.txt',wl);
    m = dlmread(fn);
    plot(m(:,1),m(:,2))
    title(sprintf('%d nm',wl))
    axis([100 1000 0 5000])
    
    k = k + 1;
end

subplot(6,7,42)
hold on
m1 = dlmread('background1.txt');
m2 = dlmread('background1.txt');
plot(m1(:,1),m1(:,2),'or')
plot(m2(:,1),m2(:,2),'xb')
legend('1','2')
title('Background')
axis([100 1000 0 5000])

saveas(gcf,'finding.png')
close all