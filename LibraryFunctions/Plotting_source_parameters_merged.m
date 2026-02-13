%% Plot slip

f = figure('Color', 'w','Units', 'inches', 'Position', [1, 1, 11, 6]);  

subplot(2,2,3)
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,slip);
axis normal;
cb=colorbar;
cb.Label.String = sprintf('     [m]'); 
cb.Label.Position = [0, 4.1, 0];      % Move label above colorbar
cb.Label.Rotation = 0;                % Make label horizontal
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 12;
xlabel('L [km]'); ylabel('W [km]');
title('Slip Hybrid')
set(gca,'fontSize',12)
clim([0 4.0])

subplot(2,2,1)
% [Nz,Nx]=size(Deterministic);
% [x,z] = meshgrid((0:Nx-1) * 1000, (0:Nz-1) * 625); 
% imagesc(x(:)./1000, z(:)./1000,Deterministic);
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,Deterministic_resampled);
axis normal;
cb=colorbar;
cb.Label.String = sprintf('     [m]'); 
cb.Label.Position = [0, 5.6, 0];      % Move label above colorbar
cb.Label.Rotation = 0;                % Make label horizontal
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 12;
xlabel('L [km]'); ylabel('W [km]');
title('Slip Deterministic')
set(gca,'fontSize',12)
clim([0 5])

subplot(2,2,[2 4])
 loglog(radial_spectra.k,radial_spectra.Pkrad_S,'k-',radial_spectra.k,radial_spectra.Pkrad_D,'b-',radial_spectra.k,radial_spectra.Pkrad_T,'g--','linewidth',1.5);
 hold on; plot([2*pi/wavelength_cut 2*pi/wavelength_cut],[max(radial_spectra.Pkrad_T) max(radial_spectra.Pkrad_T)/1000],'r--')
 legend('Stochastic','Deterministic','Hybrid','k_c','FontSize',12)
 grid on
 ylim([1e-1,1e4])
 xlabel('Radial wavenumber [rad/m]'); ylabel('Spectral amplitude')
  set(gca,'FontSize',12)

save_fig_name = sprintf('%s%s%s.png','Slip_',rup.Filename);
% exportgraphics(fig, [cd,'\',rup.Filename,'\',save_fig_name,'.pdf'], 'ContentType', 'vector');
saveas(gcf,[cd,'\',rup.Filename,'\',save_fig_name],'png');

%% Plot PDF for hypocenter

% fig=figure('Units', 'inches', 'Position', [1, 1, 7, 4.1]);
% imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,slip_pos);
% axis normal;
% cc=colorbar;
% ylabel(cc,'Hypocenter PDF','FontSize',12);
% xlabel('L [km]'); ylabel('W [km]');
% set(gca,'fontSize',12)
% hold on; plot(rup.shyp/1000,rup.dhyp/1000,'marker','pentagram','color','w','MarkerSize',8,'LineWidth',2)
% 
% save_fig_name = sprintf('%s%s%s.png','Hypo_',rup.Filename);
% saveas(gcf,[cd,'\',rup.Filename,'\',save_fig_name],'png');


%% Plot kinematic rupture parameters

f = figure('Color', 'w','Units', 'inches', 'Position', [1, 1, 10, 7]);  

t = tiledlayout( 2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax=nexttile(t);
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,rup.slip.dist{1,1});
axis normal;
cb = colorbar(ax);
% clim([0 2.6])
cb.Label.String = sprintf('    [m]'); 
cb.Label.Position = [0, max(rup.slip.dist{1,1}(:))+0.1, 0];      % Move label above colorbar
cb.Label.Rotation = 0;                % Make label horizontal
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 12;
% xlabel('L [km]'); 
ax.XTickLabel = [];
ylabel('W [km]');
title('a)                               Slip')
ax.TitleHorizontalAlignment = 'left'; 
% xlabel('L [km]');
% ylabel('W [km]');
set(ax,'fontSize',12)

ax=nexttile(t);
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,rup.psv2.dist{1,1});
axis normal;
cb = colorbar(ax);
cb.Label.String = sprintf('    [m/s]'); 
cb.Label.Position = [0, max(rup.psv2.dist{1,1}(:))+0.15, 0];      % Move label above colorbar
cb.Label.Rotation = 0;                % Make label horizontal
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 12;
ax.XTickLabel = [];
ax.YTickLabel = [];
title('b)                             V_{max}')
ax.TitleHorizontalAlignment = 'left'; 

set(gca,'fontSize',12)

ax=nexttile(t);
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,rup.Vr.dist{1,1});
axis normal;
clim([0.3 0.95]);
cb = colorbar(ax);
cb.Label.String = sprintf('    [-]'); 
cb.Label.Position = [0, 0.968, 0];      % Move label above colorbar
cb.Label.Rotation = 0;                % Make label horizontal
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 12;
xlabel('L [km]'); 
ylabel('W [km]');
% ax.XTickLabel = [];
% ax.YTickLabel = [];
title('c)                          V_{rup} / V_{S}')
ax.TitleHorizontalAlignment = 'left'; 
set(ax,'fontSize',12)
hold on
fixedLevels = 0:1:20;
[C,h]=contour(rup.x_loc_matrix./1000, rup.z_loc_matrix./1000,rup.rupT.dist{1,1}, fixedLevels, 'k');
clabel(C, h, 'FontSize', 8, 'Color', 'k');
plot(rup.shyp/1000,rup.dhyp/1000,'kp','MarkerSize',8,'MarkerFaceColor','black')

ax=nexttile(t);
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,rup.risT.dist{1,1});
axis normal;
cb = colorbar(ax);
clim([min(rup.risT.dist{1,1}(:)), 7]);
cb.Label.Position = [0, 7.15, 0];      % Move label above colorbar
cb.Label.Rotation = 0;                % Make label horizontal
cb.Label.VerticalAlignment = 'bottom';
cb.Label.String = sprintf('    [s]'); 
cb.Label.FontSize = 12;
xlabel('L [km]'); 
ax.YTickLabel = [];
title('d)                             T_{rise}')
ax.TitleHorizontalAlignment = 'left'; 
set(ax,'fontSize',12)

save_fig_name = sprintf('%s%s%s.png','Rupture_',rup.Filename);
saveas(gcf,[folder_name,'\',save_fig_name],'png');

%% Plot for peak time

% figure('Units', 'inches', 'Position', [1, 1, 8, 4]);
% t = tiledlayout( 1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% ax=nexttile;
% imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,rup.psv.dist{1,1});
% axis normal;
% cb = colorbar(ax);
% cb.Label.String = sprintf('V_{max} [m/s]'); 
% cb.Label.FontSize = 12;
% % ax.XTickLabel = [];
% % ax.YTickLabel = [];
% xlabel('L [km]'); 
% ylabel('W [km]');
% set(gca,'fontSize',12)
% 
% ax=nexttile;
% imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,rup.PT.dist{1,1});
% axis normal;
% cb = colorbar(ax);
% cb.Label.String = sprintf('Peak time [-]'); 
% cb.Label.FontSize = 12;
% % ax.XTickLabel = [];
% ax.YTickLabel = [];
% xlabel('L [km]'); 
% % ylabel('W [km]');
% set(gca,'fontSize',12)
% 
% save_fig_name = sprintf('%s%s%s.png','PeakTime',rup.Filename);
% saveas(gcf,[folder_name,'\',save_fig_name],'png');

%% Plot of roughness and integrated stress drop
figure
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t)
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,slip);
axis equal tight;
cc=colorbar;
ylabel(cc,'Slip [m]','FontSize',12);
ax = gca;
ax.XTickLabel = [];
% xlabel('L [km]'); 
ylabel('W [km]');
set(ax,'fontSize',12)

nexttile(t)
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,ISF_norm);
axis equal tight;
cc=colorbar;
ylabel(cc,'ISD [-]','FontSize',12);
ax = gca;
ax.XTickLabel = [];
% xlabel('L [km]'); 
ylabel('W [km]');
set(ax,'fontSize',12)

nexttile(t)
imagesc(rup.x_loc_matrix(:)./1000, rup.z_loc_matrix(:)./1000,Roughness);
axis equal tight;
cc=colorbar;
ylabel(cc,'Roughness [m]','FontSize',12);
xlabel('L [km]'); ylabel('W [km]');
set(gca,'fontSize',12)

save_fig_name = sprintf('%s%s%s.png','Roughnes_',rup.Filename);
saveas(gcf,[folder_name,'\',save_fig_name],'png');

%% Spectrum

figure('Units', 'inches', 'Position', [1, 1, 7, 4]);
subplot(1,2,1)
plot(rup.mrf.time{1,1},rup.mrf.mrf{1,1},'k-','LineWidth',1.5);
xlabel('t [s]','fontSize',12); 
ylabel('Moment rate [N-m]','fontSize',12);
xlim([0 max(rup.rupT.dist{1,1}(:))+3])
grid on
set(gca,'fontSize',12)

f=logspace(-2,1,200);
subplot(1,2,2)
%loglog(f,sum(rup.Mo_vec)./(1+(f/fc).^2),'r-','LineWidth',1.5); hold on
loglog(rup.mrf.freq{1,1},rup.mrf.fas{1,1},'k--','LineWidth',1.2);
xlabel('Frequency [Hz]','fontSize',12); 
ylabel('FAS [N-m/s/Hz]','fontSize',12);
xlim([0.01 10])
ylim([1e14 1e20])
grid on
set(gca,'fontSize',12)

save_fig_name = sprintf('%s%s%s.png','MomentRate_',rup.Filename);
saveas(gcf,[folder_name,'\',save_fig_name],'png');