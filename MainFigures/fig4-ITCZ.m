%(5N-5S, 170W-120W)
ts_d = ncread('Regridded_tas_18models.nc','hose');
ts_ctrl = ncread('Regridded_tas_18models.nc','control');
pr_d = ncread('Regridded_precip_18models.nc','hose');
pr_ctrl = ncread('Regridded_precip_18models.nc','control');

ts_hose = ts_ctrl+ts_d;
pr_hose = pr_ctrl+pr_d;
lat = ncread('Regridded_tas_18models.nc','lat');
lon = ncread('Regridded_tas_18models.nc','lon');
lon(lon<0) = lon(lon<0)+360;
lon=lon([73:end 1:72]);

ts_hose=ts_hose([73:end 1:72],:,:);
ts_ctrl=ts_ctrl([73:end 1:72],:,:);
ts_d=ts_d([73:end 1:72],:,:);
pr_hose=pr_hose([73:end 1:72],:,:);
pr_ctrl=pr_ctrl([73:end 1:72],:,:);
pr_d=pr_d([73:end 1:72],:,:);

%
listID = fopen(['Regridded_model_names_18models.txt']);
list = textscan(listID,'%s', 'Delimiter','\n'); % A is the entire thing
fclose(listID);
%
for i = 1: length(list{1})
    models{i,1} = list{1}{i};
end

%%

ind_lat_trop=find(lat>-20 & lat<20);

ind_lat_na1 = interp1(lat, 1:length(lat),10,'nearest','extrap');
ind_lat_na2 = interp1(lat, 1:length(lat),20 ,'nearest','extrap');
ind_lon_na1 = interp1(lon, 1:length(lon), 360-85,'nearest','extrap');
ind_lon_na2 = interp1(lon, 1:length(lon), 345,'nearest','extrap');


%% scatterplot
for i = 1:length(models)
    itcz_ctrl(:,i)=(lat(ind_lat_trop)'*squeeze(pr_ctrl(:,ind_lat_trop,i))'./squeeze(sum(pr_ctrl(:,ind_lat_trop,i),2))')';
    itcz_hose(:,i)=(lat(ind_lat_trop)'*squeeze(pr_hose(:,ind_lat_trop,i))'./squeeze(sum(pr_hose(:,ind_lat_trop,i),2))')';
end
itcz_shift = itcz_hose-itcz_ctrl;
itcz_shift_na = mean(itcz_shift(lon<15 | lon>320,:),1);

ts_d_tna = squeeze(mean(mean(ts_d(lon>300 & lon < 340, lat>0 & lat<30,:),2),1)); %tropical north atlantic

%% Plot colored ITCZs

% colormap
ice = cmocean('ice',15);
range = -4.5+0.15:0.3:0-0.15;
for i = 1: length(models)
   % ind_color(i) = interp1(x, 1:length(x),(ts_d_tna(i)-tmin)/(tmax-tmin),'nearest','extrap');
    ice_color(i) = interp1(range, 1:15, ts_d_tna(i),'nearest','extrap');
end
%%
load topo
lon_topo = 0.5:359.5;
lat_topo = -89.5:89.5;
%%
figure(9)
clf
hold on
for i = 1:length(models)
    cc = (ts_d_tna(i)+2)/(4)*15;
    %ind_cc = interp1(1:15, 1:length(1:15), cc,'nearest','extrap');
    
    plot(lon,itcz_shift(:,i),'color',ice(ice_color(i),:),'linewi',1)
end
colormap(ice)
%colormap(jet([16]))
%caxis([tmin tmax])
c = colorbar;
c.Label.String = 'Changes in Northern Tropical Atlantic temperature';
c.Label.FontSize = 12;
c.Color = 'k';
c.Ticks = -4.5:0.9:0;
caxis([-4.5 0])
plot([0 360], [0 0],'k-','linewi',1)
xlim([0 360])
xlabel('Longitude')
ylim([-10 4])
ylabel('Meridional displacement of the ITCZ')
yyaxis right
set(gca,'FontSize',12)
plot(lon_topo, mean(topo(lat_topo<5&lat_topo>-5,:),1),'k-','linewi',1)
ylim([0 20000])
yticks([])
hold off
text(200,1000,'Pacific','HorizontalAlignment','center','FontSize',12)
text(340,1000,'Atlantic','HorizontalAlignment','center','FontSize',12)
text(70,1000,'Indian','HorizontalAlignment','center','FontSize',12)
set(gca,'FontSize',12,'XTick',[0:60:360],'TickDir','out')

%%
set(0, 'DefaultFigureWindowStyle', 'docked')

sig=0.05;
for i=1:length(lon)
    for j=1:length(lat)
        a=[itcz_shift_na',  squeeze(ts_d(i,j,:))];
        %a=[nino34_hose - nino34_ctrl,  squeeze(ts_d(i,j,:))];
        [r2,p2]=corrcoef(a,'Rows','pairwise');
        r(i,j)=r2(1,2);
        p(i,j)=p2(1,2);
    end
end
%
lon30 = lon;
lon30(lon30<=30)=lon30(lon30<=30)+360;
figure(21)
clf
%m_proj('Miller','lon',[240 390],'lat',[-60 70])
m_proj('Robinson','lon',[30 390],'lat',[-80 80])
m_pcolor(lon30([14:end,1:13]),lat,r([14:end,1:13],:)')
shading interp
%colormap(brewermap([12],'BrBg'))
colormap(flipud(brewermap([20],'RdBu')))
caxis([-1 1])
c = colorbar;
c.FontSize = 12;
c.Color = 'k';
c.Label.String = 'Correlation coefficient r';
hold on
for i=1:length(lon)
    for j=1:length(lat)
        if p(i,j)<sig
            %m_line(lon30(i),lat(j),'Marker','.','Color',[.5,.5,.5],'MarkerSize',5)
            m_line(lon(i),lat(j),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',1)

        end
    end
end
% [wlat,wlon] = borders('Continental US'); 
%     for ii=1:length(wlon)
%         m_line(wlon{ii}+360,wlat{ii},'linestyle','-','Color','k','linewi',0.1);
%     end
m_grid('fontsize',12,'linewidth',1,'ytick',-60:30:60,'xtick',00:60:360+60);
m_coast('color','k','linewidth',1);
%title('Correlation between Atlantic ITCZ shift and \DeltaT')



