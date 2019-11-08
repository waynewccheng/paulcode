close all; 
clearvars;

% Metric
U = @(p) 100*(mean(p)-std(p))/(mean(p)+std(p));
CR = @(img) (max(max(img)) - min(min(img)))/(max(max(img)) + min(min(img)));
AD = @(img) std2(img)/mean2(img);

% Load
p_rdata = 'C:\Users\Paul_Lemaillet\Documents\WholeSlideImaging\ProgMatlab\Data\110819\Flatness_Of_Field\Before_FrostedGlass_NewSourceToCondenserDist';
load([p_rdata '\int_mean_array_bb']);
before_bb = int_mean_array_bb;
load([p_rdata '\int_mean_array_g']);
before_g = int_mean_array_g;

p_rdata = 'C:\Users\Paul_Lemaillet\Documents\WholeSlideImaging\ProgMatlab\Data\110819\Flatness_Of_Field\After_FrostedGlass_NewSourceToCondenserDist';
load([p_rdata '\int_mean_array_bb']);
after_bb = int_mean_array_bb;
load([p_rdata '\int_mean_array_g']);
after_g = int_mean_array_g;

% Contrast Ratio
cr_before_bb = CR(before_bb);
cr_before_g = CR(before_g);
ad_before_bb = AD(before_bb);
ad_before_g = AD(before_g);

cr_after_bb = CR(after_bb);
cr_after_g = CR(after_g);
ad_after_bb = AD(after_bb);
ad_after_g = AD(after_g);

% Images
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 2, 1);
imagesc(before_bb - min(min(before_bb)) );colorbar;
ylabel('BB');
title(['CR = ' num2str(cr_before_bb) ' AD = ' num2str(ad_before_bb)]);
subplot(2, 2, 2);
imagesc(after_bb - min(min(after_bb)) );colorbar;
title(['CR = ' num2str(cr_after_bb) ' AD = ' num2str(ad_after_bb)]);

subplot(2, 2, 3);
imagesc(before_g - min(min(before_g)) );colorbar;
ylabel('550nm')
title(['CR = ' num2str(cr_before_g) ' AD = ' num2str(ad_before_g)]);
subplot(2, 2, 4);
imagesc(after_g - min(min(after_g)) );colorbar;
title(['CR = ' num2str(cr_after_g) ' AD = ' num2str(ad_after_g)]);

figure;
y = [cr_before_bb ad_before_bb cr_after_bb ad_after_bb; cr_before_g ad_before_g cr_after_g ad_after_g];
b = bar(y,'FaceColor','flat');
for k = 1:2:size(y,2)
    b(k).CData = [0 0 1];
end

for k = 2:2:size(y,2)
    b(k).CData = [1 0 1];
end
legend('CR', 'AD')