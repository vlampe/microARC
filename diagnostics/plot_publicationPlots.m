% create figures for paper here
% code was developed in the various plot_* scripts in the /diagnostics/
% directory

 %% load model output 
close all
clear all


modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_fullRunForPaper';

%% add paths 
addpath(genpath('~/Documents/microarc/microARC model')) % location of my work with the model
addpath(genpath('~/microARC/')) % location of model


% set directories
Directories.resultsDir  = ['~/Documents/microarc/microARC model/2nd paper/mod output/' modTag '/'];
Directories.plotDir  = ['~/Documents/microarc/microARC model/2nd paper/mod output/' modTag '/plots/publicationPlots/'];


% load model output
summer = load([Directories.resultsDir  'summer_output_' modTag '.mat']);
% structure that contains model input and output for summer set-up: 
% Data: observational data used for model optimisation
% FixedParams: constants, info about depth layers, size classes, useful indeces...
% Forc: Forcing data input as time series along model trajectories
% Params: variable model parameters, equal in summer and autumn set-up
% auxVars: more model output, derived quantities or rates... (eg cell densities)
% modData: model output in the same format as Data, for convinient comparisions (but there's probably an error in the size/sizeFull substruct)
% out: model output for N (DIN), P (Phy biomass), Z (zoo biomass), OM (Detritus)
% v0: initial model conditions

% common format:
% size(out.P) = [nSizeClass nDepth nVar nDays nTrajs]  with Var = {C, N, Chl}
% size(out.OM) = [nOMType nDepth nVar nDays nTrajs] with OMType = {DOM, POM}

autumn = load([Directories.resultsDir  'autumn_output_' modTag '.mat']);

%

min(summer.out.P, [], 'all')
min(summer.out.Z, [], 'all')
min(autumn.out.P, [], 'all')

%% calc geographic distance between sampling stations

stations = unique(table(summer.Data.scalar.Label, summer.Data.scalar.Latitude, summer.Data.scalar.Longitude, 'VariableNames', ["Label", "Lat", "Lon"]), 'rows')

stations(1,:)
[arclen,az] = (distance(stations{1, 2}, stations{1, 3}, stations{2,2}, stations{2,3}))

deg2km(distance(stations{1, 2}, stations{1, 3}, stations{2,2}, stations{2,3}))


distance(78.235, 0.004, 78.333, -0.008)


dist_km = nan(height(stations), height(stations))

for s = 1:height(stations)
    dist_km(:,s)  = deg2km(distance(stations{s, 2}, stations{s, 3}, stations{:,2}, stations{:,3}))
end



%%   METHODS



%% Fig 1: plot 2x2 plot for maps, with sampling stations


watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

% map settings
land = readgeotable("landareas.shp");
proj = projcrs(3995, Authority="EPSG"); % 'stereo',  EPSG code https://epsg.io/3995

% ticks for colorbar
Ticks = [15 74 135 196 258 319]
dates = datetime(autumn.Forc.t,'ConvertFrom','datenum');
TickLabels = datestr(dates(Ticks, 1));
colormap("turbo")

t = tiledlayout(2,2, 'TileSpacing','compact')
%title(t, ['Model trajectories'])

for j = 1:length(seasons)

    season = seasons{j};
    
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
             
        
        % get coordinates
        lat = alldat.(season).Forc.y(:, wm_index);
        lon = alldat.(season).Forc.x(:, wm_index);

        % get dates and doy
        doy = yearday(alldat.(season).Forc.t(:, wm_index)); 
        dates = datetime(alldat.(season).Forc.t(:, wm_index),'ConvertFrom','datenum');

        % get sampling stations
        scalar_stations = unique(table(alldat.(season).Data.scalar.Latitude, ...
            alldat.(season).Data.scalar.Longitude, VariableNames=["lat", "lon"]));

        
        % plot maps:

        nexttile
        m = newmap(proj)

         for itraj = 1:size(lat, 2)
            
            trajtable = table(lat(:,itraj), lon(:,itraj), doy(:,itraj), VariableNames=["lat" "lon" "doy"]);
            trajgeotable = table2geotable(trajtable);
            p = geoplot(trajgeotable, ColorVariable='doy', MarkerFaceAlpha=0.5, MarkerSize=0.01);
            if itraj == 1; hold on; end
         end
         % add land
        p2 = geoplot(land(2:end, :)) % remove antarctica, else the projection produces a bug
       % p2.FaceColor = [0 0.4470 0.7410]
        p2.FaceColor = [0.5 0.5 0.5]
        % add start point of trajs
        trajtable = table([lat(1,:) lat(end,:)]', [lon(1,:) lon(end,:)]', [doy(1,:) doy(end,:)]', VariableNames=["lat" "lon" "doy"]);
        trajgeotable = table2geotable(trajtable);
        %pStart = geoplot(trajgeotable(1:height(trajtable)/2, :), '*', ColorVariable='doy', MarkerFaceAlpha=0.5, MarkerSize=3);  % oder default shape in size 10
        % add scalar sampling stations
        stationgeotable = table2geotable(scalar_stations); 
        pStations = geoplot(stationgeotable, "+", MarkerSize=3, MarkerEdgeColor="k", MarkerEdgeAlpha=0.6) %

       
        geolimits([65 85], [-20 20])
        hold off
        title([season ' setup, ' wm ' trajectories'])
        %cb = colorbar('Ticks', Ticks, 'TickLabels', TickLabels);
        



    end

end
% styling
t.Padding = 'compact'
t.TileSpacing = 'compact'

% common colorbar
cb = colorbar('Ticks', Ticks, 'TickLabels', TickLabels);
cb.Layout.Tile = "south"


colormap("turbo")

set(gcf, 'Position', [-899  -172   697   464])


% and save 
%exportgraphics(t, [Directories.plotDir 'Fig1_map_doy_seperate.png'])
% fig = gcf;
% ax = fig.CurrentAxes;
% 
% exportgraphics(ax, [Directories.plotDir 'Fig1map_doy.pdf'])
% better save by hand

%clearvars -except autumn summer Directories modTag pSetName  


%% Fig 2: plot with surface PAR and surface temp 

seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 


t = tiledlayout(2,2)
%title(t, 'Surface PAR and temperature for all trajectories')
cmap = lines(2);

for j = 1:length(seasons)

    season = seasons{j};

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];

    % PAR surf
    nexttile
    p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ...
        squeeze(alldat.(season).Forc.PARsurf(:,:,strcmp(alldat.(season).Forc.waterMass, 'Arctic'))), ...
        'Color', cmap(1,:), 'DisplayName', 'Arctic')
    
    ylabel('PAR (µEin d^{-1} m^{-2})')
    xlim([datetime('2018-01-01') datetime('2018-12-31')])
    ylim([0,5e7])
    title(['PAR ' season ' setup'])
    hold on
    p1 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ...
        squeeze(alldat.(season).Forc.PARsurf(:,:,strcmp(alldat.(season).Forc.waterMass, 'Atlantic'))), ...
        'Color', cmap(2,:), 'DisplayName', 'Atlantic')
  %  legend([p1(1) p2(1)], {'Atlantic', 'Arctic'})
    patch([data_dates fliplr(data_dates)], [min(ylim)*[1 1] max(ylim)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
    hold off
    set(gca,'children',flipud(get(gca,'children')))

    % Temp surf
    nexttile
    p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ...
        squeeze(alldat.(season).Forc.T(1,:,strcmp(alldat.(season).Forc.waterMass, 'Arctic'))), ...
        'Color', cmap(1,:), 'DisplayName', 'Arctic')
    
    ylabel('Temp (°C)')
    xlim([datetime('2018-01-01') datetime('2018-12-31')])
    ylim([-2,8.5])
    title(['Temp ' season ' setup'])
    hold on
    p1 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ...
        squeeze(alldat.(season).Forc.T(1,:,strcmp(alldat.(season).Forc.waterMass, 'Atlantic'))), ...
        'Color', cmap(2,:), 'DisplayName', 'Atlantic')
    %legend([p1(1) p2(1)], {'Atlantic', 'Arctic'})
    patch([data_dates fliplr(data_dates)], [min(ylim)*[1 1] max(ylim)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
    hold off
    set(gca,'children',flipud(get(gca,'children'))) % change order so that patch goes to background
end 

leg = legend([p1(1) p2(1)], {'Atlantic', 'Arctic'}, 'NumColumns', 2)
leg.Layout.Tile = "south"


t.TileSpacing = "compact";
t.Padding = "compact"

set(gcf, 'Position', [4   372   808   605])
saveas(t, [Directories.plotDir 'Fig2_ts_forcPARsurf_Tempsurf_watermass.png'])

clearvars -except autumn summer Directories modTag pSetName 


%% Fig 2 ALTERNATIVE: plot with surface PAR and surface temp 
%   and with secondary x and y axis

seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

fig = figure
% t = tiledlayout(2,2)  % makes secondary axes overlap
% t.TileSpacing = "compact";
% t.Padding = "compact"


% set(fig, 'Position', [4   372   808   605])
%title(t, 'Surface PAR and temperature for all trajectories')

set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 19, 16])

cmap = lines(2);

sp=0;

for j = 1:length(seasons)

    season = seasons{j};

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
    referenceDates = data_dates;

    dates = datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum');

    dates_l = datetime("2018-01-01"):datetime("2018-12-31")

    dateDiff = NaN(size(dates_l));
    for k=1:length(dateDiff)
        dat = dates_l(k);
        if dat <= min(referenceDates)
            datDist = days(dat - min(referenceDates));
        elseif dat >= max(referenceDates)
            datDist = days(dat - max(referenceDates));
        else
            datDist = 0;
        end
        dateDiff(k) = datDist;
    end
    % find every 7th date around the target time
    % dat_i = sort([find(dates == referenceDates(1)):-7:0,...
    %     find(dates == referenceDates(2)):7:length(dates)])
    %
    % dat_i = dat_i(5:end) % drop 4 first labels, they overelap with exponent from y axis
    %

    dat_i = find(day(dates_l) == 1);

    % nexttile(t)
    sp = sp+1;
    subplot(2,2,sp)
    Pos{sp} = get(gca, "Position");
    if strcmp(season, "autumn"); lower = 0.1; else; lower = Pos{sp}(2); end
    set(gca, 'Position', [0.0647, lower, 0.4, Pos{sp}(4)])

    patch([data_dates fliplr(data_dates)], [0*[1 1] 5*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.8,'edgealpha',1)
    hold on

    p1 = plot(dates, ...
        squeeze(alldat.(season).Forc.PARsurf(:,:,strcmp(alldat.(season).Forc.waterMass, 'Atlantic')))/1e7, ...
        'Color', cmap(2,:), 'DisplayName', 'Atlantic')
    ylabel('PAR (\cdot 10^7 µEin d^{-1} m^{-2})')
    ylim([0,5])

    p2 = plot(dates, ...
        squeeze(alldat.(season).Forc.PARsurf(:,:,strcmp(alldat.(season).Forc.waterMass, 'Arctic')))/1e7, ...
        'Color', cmap(1,:), 'DisplayName', 'Arctic')


    % secondary x axis
    ax1 = gca;
    set(ax1, 'FontSize', 10)
     xtickformat(ax1, "MMMMM")
    ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'YTick', []);

    xline(ax2, dates_l(dat_i), 'color', 'none');
    %xline(ax2, dates_l(dat_i), 'color', 'k');
    linkaxes([ax1, ax2]);

    % Customize the secondary x-axis labels
    xticks(ax2, dates_l(dat_i));
    xticklabels(ax2, dateDiff(dat_i))
    xlabel(ax2, 'Temporal distance from target region (d)')
    ax2.XAxis.FontSize = 8;
    % ax2.XLabel.FontSize = 9;
    % add minor labels on the 15th of each month, without label.
    set(ax1,'XMinorTick','on')
    xAx1 = get(ax1,'XAxis');
    xAx1.MinorTickValues=dates_l(dat_i)+15;
    
    xlim(ax1, [datetime('2018-01-01') datetime('2018-12-31')])
    text(10, 4.7, ['PAR ' season ' setup'], 'FontSize',10, 'FontWeight','bold')

    hold off



    % Temp surf
    % nexttile(t)
    sp = sp+1;
    subplot(2,2,sp)
    Pos{sp} = get(gca, "Position");
    % set(gca, "Position", [0.52 Pos{sp}(2:4)]);
    if strcmp(season, "autumn"); lower = 0.1; else; lower = Pos{sp}(2); end
    set(gca, 'Position', [Pos{sp}(1), lower, 0.4, Pos{sp}(4)])


    patch([data_dates fliplr(data_dates)], [-2*[1 1] 10*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.8,'edgealpha',1)
    hold on
    p1 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ...
        squeeze(alldat.(season).Forc.T(1,:,strcmp(alldat.(season).Forc.waterMass, 'Atlantic'))), ...
        'Color', cmap(2,:), 'DisplayName', 'Atlantic')
    ylabel('Temp (°C)')

    ylim([-2,8.5])


    p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ...
        squeeze(alldat.(season).Forc.T(1,:,strcmp(alldat.(season).Forc.waterMass, 'Arctic'))), ...
        'Color', cmap(1,:), 'DisplayName', 'Arctic')


    % secondary x-axis
    ax1 = gca;
    set(ax1, 'FontSize', 10)
    xtickformat(ax1, "MMMMM")
    ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'YTick', []);

    xline(ax2, dates_l(dat_i), 'color', 'none');
   % xline(ax2, dates_l(dat_i), 'color', 'k');
    linkaxes([ax1, ax2]);
    % Customize the secondary x-axis labels
    xticks(ax2, dates_l(dat_i));
    xticklabels(ax2, dateDiff(dat_i))
    xlabel(ax2, 'Temporal distance from target region (d)')
    ax2.XAxis.FontSize = 8;
    % ax2.XLabel.FontSize = 9;
    % add minor labels on the 15th of each month, without label.
    set(ax1,'XMinorTick','on')
    xAx1 = get(ax1,'XAxis');
    xAx1.MinorTickValues=dates_l(dat_i)+15;

    pause(0.5)
    xlim(ax1, [datetime('2018-01-01') datetime('2018-12-31')])
    text(10, 7.8, ['Temp ' season ' setup'], 'FontSize',10, 'FontWeight','bold')
    hold off
   
end 

leg = legend([p1(1) p2(1)], {'Atlantic', 'Arctic'}, 'NumColumns', 2)
leg.Position = [0.405 0.005 0.1683 0.0248] % set position manually



% saveas(gcf, [Directories.plotDir 'Fig2_ts_forcPARsurf_Tempsurf_watermass.png'])

savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig2_ts_forcPARsurf_Tempsurf_watermass_rev1.png']
saveas(fig, savepath)


clearvars -except autumn summer Directories modTag pSetName 

%% fig 3: DIN WOA and model

% WOA DIN profiles

% plot trajectory means of WOA data extracted along the model trajectories
%  for model validation


% _________________________________________________________________________

% 4 x 2 tiled layout
% summer arctic model   | summer arctic WOA
% summer atlantic model | summer atlantic WOA
% autumn arctic model   | autumn arctic WOA
% autumn atlantic model | autumn atlantic WOA
% 

% load WOA data
WOAtable.summer = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/monthlyTrajNSummer.csv'); % mu mol kg^-1
WOAtable.autumn = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/monthlyTrajNAutumn.csv');

% Seawater Potential Density: 1025 kg/m^3  see https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
WOAtable.summer.mean_N = WOAtable.summer.mean_N * 1025 / 1000; % convert mumol/kg to mmol/m3
WOAtable.autumn.mean_N = WOAtable.autumn.mean_N * 1025 / 1000;

% reassign model output
alldat.summer = summer;
alldat.autumn = autumn; 


seasons = [{'summer'}, {'autumn'}];
watermasses = [{'Arctic'}, {'Atlantic'}];


clabelMod = 'DIN (mmol N m^{-3})';
% clabelObs = 'DIN (µmol N kg^{-1})';
clabelObs = 'DIN (mmol N m^{-3})';
ul = 20;
ll = 0.2;

TickLabels = [0.2 1 2 4 8 16];
Ticks = log10(TickLabels);

x_ticks = datetime(['2018-01-01'; '2018-02-01'; '2018-03-01';...
            '2018-04-01'; '2018-05-01'; '2018-06-01';...
            '2018-07-01'; '2018-08-01'; '2018-09-01';...
            '2018-10-01'; '2018-11-01'; '2018-12-01']) 


fig = figure

set(fig, 'Units', 'centimeters', 'Position', [0, 0, 18, 21.5])
fontsize(fig, 10, "points")

t = tiledlayout(4,2)
t.TileSpacing = 'compact';
t.Padding = 'compact';
%title(t, 'DIN concentration, Model vs World Ocean Atlas')

for j = 1:length(seasons)

    season = seasons{j};
    
    % restructure WOAtable to array: [depth time traj]
    WOA = WOAtable.(season);
    Traj = unique(WOA.Traj);
    Month = unique(WOA.month);
    Depth = unique(WOA.depth);

    nTraj = length(Traj);
    nMonth = length(Month);
    nDepth = length(Depth);

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
    
    WOA = reshape(WOA.mean_N, [nDepth nMonth nTraj]);
    
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        wmindex = strcmp(alldat.(season).Forc.waterMass, wm);
        
        % ____ MODEL____
        % pick Model Output 
        data = alldat.(season).out.N(:, :, :, wmindex);
        % average over trajectories
        data = squeeze(mean(data, 4, 'omitnan')); 
        %size(data)

        data10 = data;
        data10(10, :) = data(9,:);  % append copy of last row so pcolor can use depth edges on y axis
        test = data(1,:);
        test(2:11,:) = data10;

        % plot 
        nexttile
        %p = pcolor(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), alldat.(season).FixedParams.zw, log10(data10))
        % set(p, 'EdgeColor', 'none'); 
        % p = contourf(alldat.(season).Forc.t(:,1), alldat.(season).FixedParams.zw, ...
        %     data10,'ShowText','on', 'LevelStep',2)
         
        p = contourf(alldat.(season).Forc.t(:,1), [0; alldat.(season).FixedParams.z; -300], ...
            test,'-w', 'ShowText','on', 'LevelStep',1, 'TextStep', 2, 'EdgeAlpha', 0.3)
        
        if i+j == 4; 
            %c = colorbar('Location', 'southoutside', 'Ticks', Ticks, 'TickLabels', TickLabels); 
            c = colorbar('Location', 'southoutside', 'Ticks', TickLabels); 
            c.Label.String = clabelMod;
           
        end;  
        %caxis(log10([ll ul]))
        caxis([ll ul])
        set(gca,'Layer', 'top', 'ColorScale','log');
      %  xtickformat('MMM')
        
      xticks(datenum(x_ticks))
      datetick('x', 'm', 'keepticks', 'keeplimits')
        

        ylabel('Depth (m)')
        hold on 
       % contour(yearday(alldat.(season).Forc.t(:,1)), alldat.(season).FixedParams.z, data, ...
       %     'k', 'ShowText','on', 'LevelStep', 2)
        patch(datenum([data_dates fliplr(data_dates)]), [min(ylim)*[1 1] max(ylim)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.7)
        hold off
        xlim(datenum([datetime('2018-01-01') datetime('2018-12-31')]))
        ylim([-200, 0])
        %
        title([season ' setup, ' wm ': MODEL'])
        
        
        %_____ WOA ______
        
        WOAmeanWM = mean(WOA(:,:,wmindex), 3); % get mean of all trajs of this watermass
        

        nexttile
        % contourf(Month, -Depth(1:29), WOAmeanWM(1:29,:),'ShowText','on')
        p = contourf(datenum(datetime(2018, Month, 15)), -Depth(1:29), WOAmeanWM(1:29,:), ...
            '-w', 'ShowText','on', 'LevelStep',1, 'TextStep', 2, 'EdgeAlpha', 0.3)
        hold on
        patch(datenum([data_dates fliplr(data_dates)]), [min(ylim)*[1 1] max(ylim)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.7)
        hold off
        if i+j == 4; 
            c = colorbar('Location', 'southoutside', 'Ticks', TickLabels); 
            c.Label.String = clabelObs;
           % c.Ruler.TickLabelFormat = '%.2f';
        end; 
        
        caxis([ll ul])
        set(gca,'ColorScale','log');
        ylabel('Depth (m)')
        xlim(datenum([datetime('2018-01-01') datetime('2018-12-31')]))
      
        ylim([-200, 0])
        xticks(datenum(x_ticks))
        datetick('x', 'm', 'keepticks', 'keeplimits')
        

        title([season ' setup, ' wm ': WOA'])

    end

end
pause(0.5)
colormap(jet())


% set(gcf, 'Position', [1441 471 791 950]);

% saveas(t, [Directories.plotDir 'Fig3_modWOA_DIN.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig3_modWOA_DIN_rev1.png']
saveas(fig, savepath)

% close all
% clearvars -except autumn summer Directories modTag pSetName

%% Fig 4: Chl a remote sensing and model

% SAT Chl trajectories

% load RS data
rsChl.summer = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/monthlyTrajChlSummer.csv'); % mu mol kg^-1
rsChl.autumn = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/monthlyTrajChlAutumn.csv');

% reassign model output
alldat.summer = summer;
alldat.autumn = autumn; 

seasons = [{'summer'}, {'autumn'}];
watermasses = [{'Arctic'}, {'Atlantic'}];


ylab = 'Phy Chl a (mg m^{-3})';
cmap = lines(2);

depthLayer = [1:5];
fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 16])
fontsize(fig, 10, 'points')

t = tiledlayout(2,2)
t.TileSpacing = 'compact';
t.Padding = 'compact';
%title(t, 'Surface Chl a concentration, Model vs Remote sensing (upper 46 m)')

for j = 1:length(seasons)

    season = seasons{j};
    
    % restructure WOAtable to array: [depth time traj]
    rs = rsChl.(season);
    Traj = unique(rs.Traj);
    Month = unique(rs.month);
    % Depth = unique(rs.depth);

    nTraj = length(Traj);
    nMonth = length(Month);
    % nDepth = length(Depth); gibt nur 1 
    zwidth = alldat.(season).FixedParams.zwidth(depthLayer); 
    
    rs = reshape(rs.mean_chl, [1 nMonth nTraj]);
    
    ppChlindex = alldat.(season).FixedParams.PP_Chl_index;

    modDates = datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum');
    modMonths = month(modDates); 

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        wmindex = strcmp(alldat.(season).Forc.waterMass, wm);
        
        % ____ MODEL____
        % pick Model Output 
        
        mod = alldat.(season).out.P(:,depthLayer,  ppChlindex, :, wmindex); % mg Chl a mm^-3
        % % average over trajectories
        % mod = squeeze(mean(mod, 5, 'omitnan')); 
        
        % average over depth
        %mod = squeeze(mean(mod, 2, 'omitnan')); %% this is an error, it should be weighed by zwidth
        mod = sum((mod .* zwidth'./sum(zwidth)), 2); 
        % sum over size classes
        mod = squeeze(sum(mod, 1));
        
        % average over trajs
        mod_mean = mean(mod, 2);
        
        %size(mod)
          % get monthly means of mod
        T = table;
        T.mod = mod_mean;
        T.modMonths = modMonths;
        G = findgroups(T.modMonths);
        Ta = groupsummary(T, 'modMonths', 'mean');
        Ta.mean_mod
        
        
        % ____ remote sensing ____
        % pick rs Output
        
        rs_wm = rs(:,:,wmindex);
        % average over Trajectories
        rs_wm = mean(rs_wm, 3, 'omitnan');
       
        
        % plot 
        nexttile

        if strcmp(wm, "Arctic"); line_opacity = 0.5; else; line_opacity = 0.3; end

        % line_opacity = 0.2

        patch([data_dates fliplr(data_dates)], [0 0 7 7], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.9,'edgealpha',0.1)
        hold on

        p = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod, 'Color', [cmap(1,:) line_opacity], 'LineWidth', 0.5, 'LineStyle', '-') % group of trajectories
        ylabel(ylab)
        xtickformat('MMMMM')

        plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod_mean, 'Color', [1 1 1], 'LineWidth', 3) % plot mean traj in white (border)
        p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod_mean, 'Color', cmap(1,:), 'LineWidth', 1.5, 'LineStyle', '-') % plot mean traj
         % plot monthly mod means
        p3 = plot(unique(datetime(2018, modMonths, 15)), Ta.mean_mod, 'o', 'Color', cmap(1,:))
         % plot rs
        p4 = plot(datetime(2018, Month, 15), rs_wm, 'o', 'Color', cmap(2,:))
        % add highlight for sampling period
        
        hold off
        title([season ' setup, ' wm])
        % legend([p2, p3, p4], {'mod.','mod. monthly avg', 'remote sensing'})
        % set(gca,'children',flipud(get(gca,'children')))
        ax = gca
        ax.Box = "on"

    end

end

linkaxes
pause(0.5)

leg = legend([p2, p3, p4], {'mod.','mod. monthly avg', 'remote sensing'})
leg.Layout.Tile = 'south';
leg.Orientation = "horizontal"; 
% set(gcf, 'Position', [1441 899 592 522]);

% saveas(t, [Directories.plotDir 'Fig4_fit_modRS_Chl_watermass_46m.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig4_fit_modRS_Chl_watermass_46m_rev1.png']
saveas(fig, savepath)

clearvars -except autumn summer Directories modTag pSetName
close all

%% Fig 4 Alternative: Chl a remote sensing and model, 8-day remote sensing data

% SAT Chl trajectories

% load RS data
rsChl.summer = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/monthlyTrajChlSummer.csv'); % mu mol kg^-1
rsChl.autumn = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/monthlyTrajChlAutumn.csv');

rsChl8d.summer = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/8dayTrajChlSummer.csv'); % mu mol kg^-1
rsChl8d.autumn = readtable('~/Documents/microarc/microARC model/Sat&WOA trajectories/8dayTrajChlAutumn.csv');


% reassign model output
alldat.summer = summer;
alldat.autumn = autumn; 

seasons = [{'summer'}, {'autumn'}];
watermasses = [{'Arctic'}, {'Atlantic'}];


ylab = 'Phy Chl a (mg m^{-3})';
cmap = lines(2);
cmap = [ 0    0.4470    0.7410; 0.8 0 0]

depthLayer = [1:5];
fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 16])
fontsize(fig, 10, 'points')

t = tiledlayout(2,2)
t.TileSpacing = 'compact';
t.Padding = 'compact';
%title(t, 'Surface Chl a concentration, Model vs Remote sensing (upper 46 m)')

for j = 1:length(seasons)

    season = seasons{j};
    
    % restructure Chl month table to array: [depth time traj]
    rs = rsChl.(season);
    Traj = unique(rs.Traj);
    Month = unique(rs.month);
    % Depth = unique(rs.depth);

    nTraj = length(Traj);
    nMonth = length(Month);
    % nDepth = length(Depth); gibt nur 1 
    zwidth = alldat.(season).FixedParams.zwidth(depthLayer); 
    
    rs = reshape(rs.mean_chl, [1 nMonth nTraj]);


    % restructure Chl 8 day table to array: [depth time traj]
    rs8d = rsChl8d.(season);
    Traj = unique(rs8d.Traj);
    Week = unique(rs8d.week);
    refDates = unique(rs8d.ref_date)+days(4.5);
    % Depth = unique(rs.depth);

    nTraj = length(Traj);
    nWeek = length(Week);
    
    rs8d = reshape(rs8d.mean_chl, [1 nWeek nTraj]);
    
    ppChlindex = alldat.(season).FixedParams.PP_Chl_index;

    modDates = datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum');
    modMonths = month(modDates); 

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        wmindex = strcmp(alldat.(season).Forc.waterMass, wm);
        
        % ____ MODEL____
        % pick Model Output 
        
        mod = alldat.(season).out.P(:,depthLayer,  ppChlindex, :, wmindex); % mg Chl a mm^-3
        % % average over trajectories
        % mod = squeeze(mean(mod, 5, 'omitnan')); 
        
        % average over depth
        %mod = squeeze(mean(mod, 2, 'omitnan')); %% this is an error, it should be weighed by zwidth
        mod = sum((mod .* zwidth'./sum(zwidth)), 2); 
        % sum over size classes
        mod = squeeze(sum(mod, 1));
        
        % average over trajs
        mod_mean = mean(mod, 2);
        
        %size(mod)
          % get monthly means of mod
        T = table;
        T.mod = mod_mean;
        T.modMonths = modMonths;
        G = findgroups(T.modMonths);
        Ta = groupsummary(T, 'modMonths', 'mean');
        Ta.mean_mod
        
        
        % ____ remote sensing ____
        % pick rs Output
        
        rs_wm = rs(:,:,wmindex);
        % average over Trajectories
        rs_wm = mean(rs_wm, 3, 'omitnan');

        rs_wm8d = rs8d(:,:,wmindex);
        % average over Trajectories
        rs_wm8d = mean(rs_wm8d, 3, 'omitnan');
       
        
        % plot 
        nexttile

        if strcmp(wm, "Arctic"); line_opacity = 0.5; else; line_opacity = 0.3; end

        % line_opacity = 0.2

        patch([data_dates fliplr(data_dates)], [0 0 7 7], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.9,'edgealpha',0.1)
        hold on

        p = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod, 'Color', [cmap(1,:) line_opacity], 'LineWidth', 0.5, 'LineStyle', '-') % group of trajectories
        ylabel(ylab)
        xtickformat('MMMMM')

        plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod_mean, 'Color', [1 1 1], 'LineWidth', 3) % plot mean traj in white (border)
        p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod_mean, 'Color', cmap(1,:), 'LineWidth', 1.5, 'LineStyle', '-') % plot mean traj
         % plot monthly mod means
        p3 = plot(unique(datetime(2018, modMonths, 15)), Ta.mean_mod, 'o', 'Color', cmap(1,:))
         % plot rs
        p4 = plot(datetime(2018, Month, 15), rs_wm, 'o', 'Color', cmap(2,:))
        % plot 8 day rs
        p5 = plot(refDates, rs_wm8d, '+', 'Color', cmap(2,:))
        
        hold off
        title([season ' setup, ' wm])
        % legend([p2, p3, p4], {'mod.','mod. monthly avg', 'remote sensing'})
        % set(gca,'children',flipud(get(gca,'children')))
        ax = gca
        ax.Box = "on"

    end

end

linkaxes
pause(0.5)

leg = legend([p2, p3, p4, p5], {'mod.','mod. monthly avg', 'rs monthly avg', 'rs 8-day avg'})
leg.Layout.Tile = 'south';
leg.Orientation = "horizontal"; 
% set(gcf, 'Position', [1441 899 592 522]);

% saveas(t, [Directories.plotDir 'Fig4_fit_modRS_Chl_watermass_46m.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig4_fit_modRS_Chl_watermass_46m_8day_rev1.png']
saveas(fig, savepath)

% clearvars -except autumn summer Directories modTag pSetName
% close all
%% Figure 5: RESULTS: fit of model data to observations: boxplots statevars

% copied and adjusted (for total POM) from plot_modelDataFit.m

%% boxplots of true (total) POM (phy + zoo + det), without water mass distinction
% with increased time windows around observations
% not constucted from modData, but from out

% PON, POC
plotVars = [{'DIN'}, {'Chl a'}, {'POC'}, {'PON'}];

seasons = [{'summer'}, {'autumn'}];
alldat.summer = summer;
alldat.autumn = autumn;

timeFrame = 7; % include simulations +- 7 days around sampling of observations

% plot Vaules over depth, sperated by watermasses
fig = figure; 
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 24])
fontsize(fig, 10, 'points')

t = tiledlayout(length(seasons), length(plotVars));
%title(t, 'Comparision of observations to corresponding model output')

for s = 1:length(seasons)
    season = seasons{s};
    
    % observation data
    dat = alldat.(season).Data.scalar.Value;
    datScaled = alldat.(season).Data.scalar.scaled_Value;

    depth = alldat.(season).Data.scalar.Depth;

    % pick model equivalents with a wider time frame
    samplingDates = unique(alldat.(season).Data.scalar.Yearday); 
    timeindexMOD = samplingDates(1)-timeFrame:samplingDates(end)+7; 
    
    for pv = 1:length(plotVars)

        nexttile

        plotVar = plotVars{pv};
        switch plotVar
            case 'DIN'
                varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'N');
                elemVar = 'N';
                ylab = {['DIN concentration'], ['(mmol m^{-3})']};

                xticks = [0:5:20]; 
                xminorticks = [0:1:20]; 

            case 'Chl a'
                varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'chl_a');
                elemVar = 'Chl';
                ylab = {['Chl concentration'], ['(mg m^{-3})']};

                xticks = [0:6];
                xminorticks = [0:0.5:6];


            case 'POC'
                varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'POC');
                elemVar = 'C';
                ylab = {['POC concentration'], ['(mmol m^{-3})']};

                xticks = [0:10:50];
                xminorticks = [0:5:50];
            case 'PON'
                varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'PON');
                ylab = {['PON concentration'], ['(mmol m^{-3})']};
                elemVar = 'N';

                xticks = [0:1:5];
                xminorticks = [0:0.5:5];
        end


        % pick obs data to plot (according to season, watermass, variable)
        plotdat = dat(varindexDAT);
        plotdatScaled = datScaled(varindexDAT);
        depthDAT = depth(varindexDAT);

        
        
        % pick model output 
        if strcmp(plotVar, 'POC')|strcmp(plotVar, 'PON') 
             
            % detOM
            POM_index = alldat.(season).FixedParams.POM_index;
            pomvarindexMOD = strcmp(alldat.(season).FixedParams.OM_nut, elemVar);
            detOM = squeeze(alldat.(season).out.OM(POM_index, :, pomvarindexMOD, timeindexMOD, :));
            %size(detOM)
    
            % phyOM
            phyvarindexMOD = strcmp(alldat.(season).FixedParams.PP_nut, elemVar); 
            phyOM = alldat.(season).out.P(:,:,phyvarindexMOD, timeindexMOD, :); 
            phyOM = squeeze(sum(phyOM, 1)); % sum over size classes
            %size(phyOM)
    
            % zooOM
            zoovarindexMOD = strcmp(alldat.(season).FixedParams.ZP_nut, elemVar); 
            zooOM = alldat.(season).out.Z(:,:,zoovarindexMOD, timeindexMOD, :); 
            zooOM = squeeze(sum(zooOM, 1)); % sum over size classes
            %size(zooOM)
    
            % add all OM sources
            plotmod = detOM + phyOM + zooOM;
            %size(plotmod)
        elseif strcmp(plotVar, 'DIN')
            plotmod = squeeze(alldat.(season).out.N(:,:, timeindexMOD, :)); 
        elseif strcmp(plotVar, 'Chl a')
            phyvarindexMOD = strcmp(alldat.(season).FixedParams.PP_nut, 'Chl'); 
            plotmod = alldat.(season).out.P(:,:,phyvarindexMOD, timeindexMOD, :); 
            plotmod = squeeze(sum(plotmod,1));
        end

        % repair: if plotdat is empty: skip to next iteration (plotting
        % of only mod results not possible here, no depthDAt vector
        % available for interpolation
        if isempty(plotdat)
            plotdat = [NaN]; 
            depthDAT = [NaN];

            % and reshape to long vector
            plotmodVector = reshape(plotmod, [numel(plotmod), 1]); 

            depthsMOD = abs(alldat.(season).FixedParams.z);
            depthsMODVector = repmat(depthsMOD, size(plotmod,2)*size(plotmod,3),1);
    
    
            % grouping variable for color
            colg = [repmat("data", length(plotdat),1) ; repmat("model", length(plotmodVector),1)]; 
    
            plotdepth = [depthDAT; depthsMODVector];


        else

    
            % interpolate to obs depths
            depthsMOD = abs(alldat.(season).FixedParams.z);
            udepthDAT = unique(depthDAT);
            plotmodINT = NaN([length(udepthDAT), size(plotmod,2), size(plotmod,3)]);
    
            % this looped interpolation is really slow, there must be
            % another way
            for i=1:size(plotmod,2)
                for j=1:size(plotmod,3)
                    plotmodINT(:, i, j) = interp1(depthsMOD, plotmod(:,i,j), udepthDAT);
                end
            end
    
            % and reshape to long vector
            plotmodVector = reshape(plotmodINT, [numel(plotmodINT), 1]); 
            depthsMODVector = repmat(udepthDAT, size(plotmodINT,2)*size(plotmodINT,3),1);
    
    
            % grouping variable for color
            colg = [repmat("data", length(plotdat),1) ; repmat("model", length(plotmodVector),1)]; 
    
            plotdepth = [depthDAT; depthsMODVector];
        end

        % use boxchart 
        % nexttile
        % boxchart(categorical(plotdepth), [plotdat; plotmodVector], 'GroupByColor', colg, ...
        %      'MarkerStyle', '+') %, 'Orientation', 'horizontal'
        b = boxchart(categorical(round(plotdepth)), [plotdat; plotmodVector], 'GroupByColor', colg, ...
             'MarkerStyle', '.') %, 'Orientation', 'horizontal'
        b(1).JitterOutliers = 'on';
        b(2).JitterOutliers = 'on';
        
        view(90,90) % turn it, because when using groupbycolor the orientation horizontal argument does not work
        xlabel('Depth (m)')
        ylabel(ylab)
        ylim([0 Inf])
        title([season ' setup']);

        % make better x ticks 
        ax = gca
        ax.YAxis.TickValues = xticks; 
        ax.YTickLabelRotation = 0; 
        ax.YAxis.MinorTick = "on"

    end % end var loop
end % end season loop


leg = legend
leg.Layout.Tile = "south"
leg.Orientation = "horizontal"


% pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
% width = pos(3);
% height = pos(4);

t.TileSpacing = "compact";
t.Padding = "compact"

% set(gcf, 'Position', [441 1 903 796])
% saveas(t, [Directories.plotDir 'Fig5_fit_modData_boxplot_stateVars.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig5_fit_modData_boxplot_stateVars_rev1.png']
saveas(fig, savepath)

clearvars -except autumn summer Directories modTag pSetName
close all




%%
% size spectra
% abundance and biovolume
concTypes = [{'abundance'}, {'biovolume'}];

alldat.summer = summer; 
alldat.autumn = autumn;

seasons = [{'summer'}, {'autumn'}];
waterMasses = [{'Arctic'}, {'Atlantic'}];

iDepth = [1:5]; 
summer.FixedParams.zw(iDepth+1)

ESDs = summer.FixedParams.PPdia_intervals;
dlog10ESD = diff(log10(ESDs(1:2))); 

x = [ESDs(1), repelem(ESDs(2:end-1),2)', ESDs(end)]; 
x2 = [x, flip(x)];
zscore = 1.96; % for 95% CI
 cmap = lines(2); % phy blue, zoo orange
 cmap = [ 0.1 0.9 0.1; 0.1 0.1 1]
cmap = [100, 221, 23;  144 12 63 ]./255 ;
% cmap = [0.3216 0.5765 0.2275; 0.5882 0.1961 0.1412]

% cmap = [0.3216 0.5765 0.2275; 0.5 0.5 0.5]


for ct = 1:length(concTypes)
    concType = concTypes{ct};
    
    switch concType
        case 'abundance'
            modVar = 'cellDensity';
            obsVar = 'cellDensity';
            miny = 1e3; 
            totunit = '(cells / m^3)';
            yunit = '(cells / m^3 log_{10}(ESD / 1 µm)^{-1})'; 
            ytext = 'cell conc. density';
        case 'biovolume'
            modVar = 'biovolume';
            obsVar = 'BioVolDensity';
            miny = 1e6; 
            totunit = '(µm^3 / m^3)';
            yunit = '(µm^3 / m^3 log_{10}(ESD / 1 µm)^{-1})';
            ytext = 'biovol. conc. density'; 
    end
    
    fig = figure
    set(fig, 'Units', 'centimeters')
    pos = get(fig, 'Position')
    set(fig, 'Position', [0, 0, 18, 14])
    fontsize(fig, 10, 'points')

    t = tiledlayout(2,2)
    t.TileSpacing = "compact";
    t.Padding = "compact"
  %  title(t, ['Observed and simulated size spectra of plankton ' concType])

    for s=1:length(seasons)
        season = seasons{s}; 

        switch season
            case 'summer'
               % iTime = [197:208];
                iTime = [190:215];
            case 'autumn'
               % iTime = [259:277]; 
                iTime = [252:284]; 
        end

        iPhyto = alldat.(season).FixedParams.phytoplankton;
        iZoo  = alldat.(season).FixedParams.zooplankton;


        for w=1:length(waterMasses)
            wm = waterMasses{w}; 
            iTraj = strcmp(alldat.(season).Forc.waterMass, wm);

            % get observation data for comparison
            if strcmp(season, 'summer')&strcmp(wm, 'Atlantic')
                scenario = 'S1';
            elseif strcmp(season, 'summer')&strcmp(wm, 'Arctic')
                scenario = 'S2';
            elseif strcmp(season, 'autumn')&strcmp(wm, 'Atlantic')
                scenario = 'S3';
            elseif strcmp(season, 'autumn')&strcmp(wm, 'Arctic')
                scenario = 'S4';
            end
            

            scenarioindex = strcmp(alldat.(season).Data.size.scenario, scenario);
            obsESDs = alldat.(season).Data.size.ESD(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'heterotroph')); 
            dlog10ESDobs = diff(log10(obsESDs(1:2)));

            obsPhy = alldat.(season).Data.size.(obsVar)(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'autotroph'));
            obsPhySE = alldat.(season).Data.size.(strcat(obsVar, 'SE'))(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'autotroph'));
            obsPhyConc = sum(obsPhy .* dlog10ESDobs); 


            obsZoo = alldat.(season).Data.size.(obsVar)(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'heterotroph'));
            obsZooSE = alldat.(season).Data.size.(strcat(obsVar, 'SE'))(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'heterotroph'));
            obsZooConc = sum(obsZoo .* dlog10ESDobs); 
            


            % pick model output
            modPhy = alldat.(season).auxVars.(modVar)(iPhyto, iDepth, iTime, iTraj); % abundance or biovol per size class
            modPhyConc = mean(sum(modPhy,1), 2:4);
            modPhy = modPhy ./ dlog10ESD;  % convert to abundance or biovol density
            modPhySD = std(modPhy, [], 2:4); % get standard deviation
            modPhySE = modPhySD ./ sqrt(prod(size(modPhy, 2:4)));  % get standard error
            modPhy = mean(modPhy, 2:4); % mean over depth, time and trajs (biovol dens per size class is left)
            
            modZoo = alldat.(season).auxVars.(modVar)(iZoo, iDepth, iTime, iTraj); 
            modZooConc = mean(sum(modZoo,1), 2:4);
            modZoo  = modZoo ./ dlog10ESD;
            modZooSD = std(modZoo, [], 2:4);
            modZooSE = modZooSD ./ sqrt(prod(size(modZoo, 2:4)));
            modZoo = mean(modZoo, 2:4); 

            modPhyRep = repelem(modPhy', 2);
            modPhySERep = repelem(modPhySE', 2);

            modZooRep = repelem(modZoo', 2);
            modZooSERep = repelem(modZooSE', 2);



            totalconcstr = {['total ' totunit ],...
                                ['Phy Mod: ' num2str(mean(modPhyConc, 'omitnan'), '%.3e')], ...
                                ['Phy Dat: ' num2str(mean(obsPhyConc), '%.3e')], ...
                                ['Zoo Mod: ' num2str(mean(modZooConc, 'omitnan'), '%.3e')], ...
                                ['Zoo Dat: ' num2str(mean(obsZooConc), '%.3e')]};



            nexttile
               inBetween = [modPhyRep - zscore.*modPhySERep, flip(modPhyRep + zscore.*modPhySERep)];
               inBetween(inBetween <= miny) = miny;
            fill(x2, inBetween, cmap(1, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            hold on
            plot(x, repelem(modPhy', 2), 'Color', cmap(1, :), 'LineWidth', 2); % PHY
           
            inBetween = [modZooRep - zscore.*modZooSERep, flip(modZooRep + zscore.*modZooSERep)];
            inBetween(inBetween <= miny) = miny;
            fill(x2, inBetween, cmap(2, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(x, repelem(modZoo', 2), 'Color', cmap(2, :), 'LineWidth', 2); % ZOOO
          
            inBetween = [obsPhy - zscore*obsPhySE; flip(obsPhy + zscore*obsPhySE)];
            inBetween(inBetween <= miny) = miny;
            fill([obsESDs; flip(obsESDs)], inBetween, cmap(1, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(obsESDs, obsPhy, 'Color', cmap(1, :), 'LineStyle', ':', 'LineWidth', 2);
            
                
            inBetween = [obsZoo - zscore*obsZooSE; flip(obsZoo + zscore*obsZooSE)];
            inBetween(inBetween <= miny) = miny;
            fill([obsESDs; flip(obsESDs)], inBetween, cmap(2, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(obsESDs, obsZoo, 'Color', cmap(2, :), 'LineStyle', ':', 'LineWidth', 2);
            
            ax = gca
            set(ax, 'YScale', 'log', 'XScale', 'log');
            title([season ' setup, ' wm ' trajectories']);
            ylim([miny, Inf])

            xlabel('cell size / ESD (µm)')
            ylabel({ytext; yunit})

            ax.XAxis.TickValues = [1 2 5 10 20 50 100 200];
            ax.XTickLabelRotation = 0;
            ax.YAxis.TickValues = [1e4 1e6 1e8 1e10 1e12];

            
            % switch concType
            %     case 'abundance'
            %         legend('', 'Phy', '', 'Zoo','Location','northeast')
            %     case 'biovolume'
            %         legend('', 'Phy', '', 'Zoo','Location','northwest')
            % end
            
        end
    end
    linkaxes

     switch concType
        case 'abundance'
            leg = legend('', 'Phy', '', 'Zoo','Location','northeast')
        case 'biovolume'
            leg = legend('', 'Phy', '', 'Zoo','Location','northwest')
     end

     leg.Layout.Tile = "south"
     leg.Orientation = "horizontal"
    % set(gcf, 'Position', [4   372   808   605])
    % saveas(t, [Directories.plotDir 'Fig6ab_spectra_' modVar '_fitAtSamplingSSpectra.png'])

    savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig6ab_spectra_' modVar '_fitAtSamplingSSpectra_rev1.png']
    saveas(fig, savepath)
end


% clearvars -except autumn summer Directories modTag pSetName

%% Fig 7: changes in community size structure

% % Stacked Area Plots
% plot integrated C biomass per size class against time

% phytoplankton

watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

maxdepth = 8; % 100m

phycolors = ["#cde0b0"; 
    "#95a97b";
    "#5e6b39";
    "#969637";
    "#cdda71";
    "#52933a";
    "#73d670";
    "#d8db33";
    "#7bdb3e"]

sizeClasses = summer.FixedParams.PPdia; 
sizeClassEdges = summer.FixedParams.PPdia_intervals;

sizeclasslabels = cell(size(sizeClasses)); 

for sc = 1:length(sizeclasslabels)
    sizeclasslabels{sc} = [num2str(round(sizeClassEdges(sc),2)) ' - ' num2str(round(sizeClassEdges(sc+1),2)) ' µm'];
end

fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 14])
fontsize(fig, 10, 'points')

t = tiledlayout(2,2)
t.TileSpacing = "compact";
t.Padding = "compact"
% title(t, ['Phy biomass per size class (depth-integrated, upper ' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])
c = 0
for j = 1:length(seasons)

    season = seasons{j};
    
    zWidth = alldat.(season).FixedParams.zwidth(1:maxdepth); 
    
    Phy_index = alldat.(season).FixedParams.phytoplankton;
    C_index = strcmp(alldat.(season).FixedParams.PP_nut, 'C');
    
    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
 
    for i = 1:length(watermasses)
        c = c+1;
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
        
        % all phy biomass 
        mod = alldat.(season).out.P(:,1:maxdepth, C_index, :, wm_index); % mmol m^-3
        
        %average over trajs
        mod = mean(mod, 5, 'omitnan');

        % integrate over depth
        mod = squeeze(sum((mod .* zWidth'),2, 'omitnan')); % mmol m^-2
        size(mod)

        % convert go g C m^-2
        mod = mod .* (12.011 * 1e-3);
        
        nexttile
        
        patch([data_dates fliplr(data_dates)], [0*[1 1] 10*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.8,'edgealpha',1, 'DisplayName', 'none')
        % p = area(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod')
       
        hold on
        
        p = area(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), flip(mod)')
        hold off
        
        ylabel('Phy C (g C m^{-2})')
        colororder(phycolors)
        xtickformat('MMMMM')
        title([season ' setup, ' wm ' trajectories'])
        ax = gca
        ax.Box = 'on'

    end
end
linkaxes

leg = legend([p(9) p(8) p(7) p(6) p(5) p(4) p(3) p(2) p(1)], sizeclasslabels,...
    'NumColumns',3)
leg.Layout.Tile = 'south';

% set(gcf, 'Position', [4   372   808   605])

% saveas(t, [Directories.plotDir 'Fig7_Phy_C_biomass_per_sizeclass_int' num2str(sum(zWidth)) 'm.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig7_Phy_C_biomass_per_sizeclass_int' num2str(sum(zWidth)) 'm_rev1.png']
saveas(fig, savepath)

% clearvars -except autumn summer Directories modTag pSetName





% zooplankton

watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

maxdepth = 8; % 100m

zoocolors = ["#4e3021",
    "#963224";
    "#884a4c";
    "#856435";
    "#c24363";
    "#cc7240";
    "#e38181";
    "#d0a34c";
    "#ceac95"]

sizeClasses = summer.FixedParams.ZPdia; 
sizeClassEdges = summer.FixedParams.ZPdia_intervals;

sizeclasslabels = cell(size(sizeClasses)); 

for sc = 1:length(sizeclasslabels)
    sizeclasslabels{sc} = [num2str(round(sizeClassEdges(sc),2)) ' - ' num2str(round(sizeClassEdges(sc+1),2)) ' µm'];
end


fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 14])
fontsize(fig, 10, 'points')

t = tiledlayout(2,2)
t.TileSpacing = "compact";
t.Padding = "compact"
%title(t, ['Zoo biomass per size class (depth-integrated, upper ' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])

for j = 1:length(seasons)

    season = seasons{j};
    
    zWidth = alldat.(season).FixedParams.zwidth(1:maxdepth); 
    
    Zoo_index = alldat.(season).FixedParams.zooplankton;
    C_index = strcmp(alldat.(season).FixedParams.ZP_nut, 'C');
    
    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];

    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
        
        % all phy biomass 
        mod = alldat.(season).out.Z(:,1:maxdepth, C_index, :, wm_index); % mmol m^-3
        
        %average over trajs
        mod = mean(mod, 5, 'omitnan');

        % integrate over depth
        mod = squeeze(sum((mod .* zWidth'),2, 'omitnan')); % mmol m^-2
        size(mod)

        % convert go g C m^-2
        mod = mod .* (12.011 * 1e-3);
        
        nexttile
        patch([data_dates fliplr(data_dates)], [0*[1 1] 3*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.8,'edgealpha',1)
        
        hold on
        p = area(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), flip(mod)')
        hold off
        ylabel('Zoo C (g C m^{-2})')
        colororder(zoocolors)
        xtickformat('MMMMM')
        title([season ' setup, ' wm ' trajectories'])

        ax = gca
        ax.Box = 'on'
    end
end
linkaxes
%leg = legend(flip(sizeclasslabels) ,'NumColumns',3)
leg = legend([p(9) p(8) p(7) p(6) p(5) p(4) p(3) p(2) p(1)], sizeclasslabels,...
    'NumColumns',3)
leg.Layout.Tile = 'south';

% set(gcf, 'Position', [4   372   808   605])

% saveas(t, [Directories.plotDir 'FigB2_Zoo_C_biomass_per_sizeclass_int' num2str(sum(zWidth)) 'm.png'])

savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'FigB2_Zoo_C_biomass_per_sizeclass_int' num2str(sum(zWidth)) 'm_rev1.png']
saveas(fig, savepath)

clearvars -except autumn summer Directories modTag pSetName


%% domintant size class plot, phy und zoo zusammen 
% sieht man ein zeitliches lag?

% which of the isze classes contains the max biomass at the time?
% bezogen auf C biomass


maxdepth = 5; % 46m, end of euphotic zone (from visual analysis of PhyC/PhyChl profiles)

planktontypes = [{'Phyto'}, {'Zoo'}]
watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

cmap = lines(2);

sizeClasses = summer.FixedParams.PPdia; 
sizeClassEdges = summer.FixedParams.PPdia_intervals;

sizeclasslabels = cell(size(sizeClasses)); 

for sc = 1:length(sizeclasslabels)
    sizeclasslabels{sc} = [num2str(round(sizeClassEdges(sc),2)) ' - ' num2str(round(sizeClassEdges(sc+1),2)) ' µm'];
end

t = tiledlayout(2,2)
t.TileSpacing = "compact";
t.Padding = "compact"
%title(t, ['Most prominent (biomass) plankton size class (depth-integrated, upper ' num2str(round(abs(alldat.summer.FixedParams.zw(maxdepth+1)))) ' m)'])

for j = 1:length(seasons)

    season = seasons{j};
    
    z = alldat.(season).FixedParams.z;
    zEdges = alldat.(season).FixedParams.zw;
    zWidth = alldat.(season).FixedParams.zwidth; 

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
        
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);


        if strcmp(wm, "Arctic"); line_opacity = 0.4; else; line_opacity = 0.1; end

        % initialise plots
        p = cell(length(planktontypes),1);
        nexttile
        patch([data_dates fliplr(data_dates)], [0*[1 1] 10*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
        hold on
        
        for ptype = 1:length(planktontypes)

            planktontype = planktontypes{ptype}; 
            switch planktontype
                case 'Phyto'
                    P_index = 'P';
                case 'Zoo'
                    P_index = 'Z';
            end
        
            C_index = strcmp(alldat.(season).FixedParams.PP_nut, 'C');
            % all biomass in size classes
            mod = squeeze(alldat.(season).out.(P_index)(:, :, C_index, :, wm_index));
            
            % sum over size class (1st dim)
            % mod = squeeze(sum(mod, 1, 'omitnan')); 
            % integrate over depth until maxdepth
            mod = mod(:,1:maxdepth, :, :) .* zWidth(1:maxdepth)'; % mmol C /m3 / d -> mmol C /m2 / d
            mod = squeeze(sum(mod, 2, 'omitnan'));
            
            % find wich size class has the most biomass in all trajs on all
            % days
            %[row, col] = find(mod == max(mod, [], 1));
            %mod = reshape(row, size(mod, [2,3]));
            [M, I] = max(mod, [], 1);
            mod = squeeze(I);

            % plot 
    
            plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod, 'Color', [cmap(ptype,:) line_opacity])
            p{ptype} = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(mod,2), 'Color', cmap(ptype,:), 'LineWidth', 1.5)
            
        end
        
        hold off
        xlabel('time')
        yticks([1:9])
        yticklabels(sizeclasslabels)
        ylabel('size class')
        title([season ' setup, ' wm ' trajectories'])
        linkaxes;
        leg = legend([p{1}(1), p{2}(1)], planktontypes);
        set(leg, 'Location', 'Southeast')

    end
end

set(gcf, 'Position', [4   372   808   605])

saveas(t, [Directories.plotDir 'Fig8_dominSizeClass_biomass_PhyZoo_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm_allTrajs_waterMass_ts.png'])

clearvars -except autumn summer Directories modTag pSetName

%% domintant size class plot, phy und zoo zusammen WITH BIOMASS BACKGROUND


maxdepth = 5; % 46m, end of euphotic zone (from visual analysis of PhyC/PhyChl profiles)

planktontypes = [{'Phyto'}, {'Zoo'}]
watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

cmap = lines(2);

sizeClasses = summer.FixedParams.PPdia; 
sizeClassEdges = summer.FixedParams.PPdia_intervals;

sizeclasslabels = cell(size(sizeClasses)); 

for sc = 1:length(sizeclasslabels)
    % sizeclasslabels{sc} = [num2str(round(sizeClassEdges(sc),2)) ' - ' num2str(round(sizeClassEdges(sc+1),2)) ' µm'];
    sizeclasslabels{sc} = [num2str(round(sizeClassEdges(sc),2)) ' - ' num2str(round(sizeClassEdges(sc+1),2))];
end

fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 14])
fontsize(fig, 10, 'points')


t = tiledlayout(2,2)
t.TileSpacing = "compact";
t.Padding = "compact"
%title(t, ['Most prominent (biomass) plankton size class (depth-integrated, upper ' num2str(round(abs(alldat.summer.FixedParams.zw(maxdepth+1)))) ' m)'])

for j = 1:length(seasons)

    season = seasons{j};
    
    z = alldat.(season).FixedParams.z;
    zEdges = alldat.(season).FixedParams.zw;
    zWidth = alldat.(season).FixedParams.zwidth; 

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
        
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);


        if strcmp(wm, "Arctic"); line_opacity = 0.4; else; line_opacity = 0.1; end

        % initialise plots
        p = cell(length(planktontypes),1);
        nexttile
        % patch([data_dates fliplr(data_dates)], [0*[1 1] 10*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        % 'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
        patch(datenum([data_dates fliplr(data_dates)]), [0*[1 1] 10*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
        
        hold on
        
        for ptype = 1:length(planktontypes)

            planktontype = planktontypes{ptype}; 
            switch planktontype
                case 'Phyto'
                    P_index = 'P';
                case 'Zoo'
                    P_index = 'Z';
            end
        
            C_index = strcmp(alldat.(season).FixedParams.PP_nut, 'C');
            % all biomass in size classes
            mod = squeeze(alldat.(season).out.(P_index)(:, :, C_index, :, wm_index));
            
            % sum over size class (1st dim)
            % mod = squeeze(sum(mod, 1, 'omitnan')); 
            % integrate over depth until maxdepth
            mod = mod(:,1:maxdepth, :, :) .* zWidth(1:maxdepth)'; % mmol C /m3 / d -> mmol C /m2 / d
            mod = squeeze(sum(mod, 2, 'omitnan'));

            if strcmp(planktontype, 'Phyto')
                bioMass = mean(mod,3, 'omitnan'); 
                % plot background layer with phytoplankton biomass

                im = imagesc(alldat.(season).Forc.t(:,1), [1:9], ...
                        log10(bioMass), 'AlphaData',0.75);

                clim(log10([0.1 150]));
                %im.AlphaData = 0.7;
                % colLabels = [0.1 0.2 0.5 1 2 5 10 20 50 100 150];
                % colTicks = log10(colLabels);
                % 
                % clim(log10([0.1 150]))
                % c = colorbar('Ticks', colTicks, 'TickLabels', colLabels, ...
                %     'Location', 'southoutside')
                % c.Label.String = 'Phy biomass [mmol C m^{-2}]'
                % 
            end
            
            % find wich size class has the most biomass in all trajs on all
            % days
            %[row, col] = find(mod == max(mod, [], 1));
            %mod = reshape(row, size(mod, [2,3]));
            [M, I] = max(mod, [], 1);
            mod = squeeze(I);

            % plot 
            
            % plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mod, 'Color', [cmap(ptype,:) line_opacity])
            % p{ptype} = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(mod,2), 'Color', cmap(ptype,:), 'LineWidth', 1.5)
            plot(alldat.(season).Forc.t(:,1), mod, 'Color', [cmap(ptype,:) line_opacity])
            p{ptype} = plot(alldat.(season).Forc.t(:,1), mean(mod,2), 'Color', cmap(ptype,:), 'LineWidth', 1.5)
            
        end
        
        % repeat patch so it is on top
        patch(datenum([data_dates fliplr(data_dates)]), [0*[1 1] 10*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
        'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
        
        hold off

        ax = gca
        ax.Box = "on"

        %xlabel('time')
        yticks([1:9])
        yticklabels(sizeclasslabels)
        ylabel('size class (µm)')

        datetick('x', 'm', 'keeplimits') % convert datenums to datetime labels
        ax.XTickLabelRotation = 0; 
        ax.XLim = (datenum([datetime('2018-01-01') datetime('2018-12-31')]))
      
        title([season ' setup, ' wm ' trajectories'])
        linkaxes;
        % leg = legend([p{1}(1), p{2}(1)], planktontypes);
        % set(leg, 'Location', 'Southeast')

    end
end

leg = legend([p{1}(1), p{2}(1)], planktontypes);
leg.Layout.Tile = "south"
% leg.Orientation = "horizontal"
% leg.Position = [0.4423, 0.0989, 0.2401, 0.0378]

colLabels = [0.1 0.2 0.5 1 2 5 10 20 50 100 150];
colTicks = log10(colLabels);
c = colorbar('Ticks', colTicks, 'TickLabels', colLabels)
c.Label.String = 'Phy biomass [mmol C m^{-2}]'
c.Layout.Tile = 'south'
c.Layout.TileSpan = [1 1]
c.Orientation = 'horizontal'
% c.Position = [0.1755, 0.022, 0.777, 0.0302]




% 
% cpos = c.Position 
% c.Position = [cpos(1), cpos(2)-0.1, cpos(3:4)]

% set(gcf, 'Position', [4   372   808   605])


cmp = colormap("summer") % summer cmap looks nicer i think
cmp = flipud(cmp)
colormap(cmp)


% saveas(t, [Directories.plotDir 'Fig8_dominSizeClass_biomass_PhyZoo_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm_allTrajs_waterMass_ts_WITHBIOMASS.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig8_dominSizeClass_biomass_PhyZoo_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm_allTrajs_waterMass_ts_WITHBIOMASS_rev1.png']
saveas(fig, savepath)



clearvars -except autumn summer Directories modTag pSetName




%% Production and export plots (needs to run to calc production and export)

% production
maxdepth = 8; % 100m, end of euphotic zone (from visual analysis of PhyC/PhyChl profiles)

watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

% map settings
land = readgeotable("landareas.shp");
proj = projcrs(3995, Authority="EPSG"); % 'stereo',  EPSG code https://epsg.io/3995


mycolormap = [165 0 38;
    215 48 39;
    244 109 67;
    253 174 97;
    254 224 144;
    224 243 248;
    171 217 233;
    116 173 209;
    69 117 180;
    49 54 149]./255;

t = tiledlayout(2,2)
title(t, ['Net Primary Production (depth-integrated, upper ' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])

for j = 1:length(seasons)

    season = seasons{j};
    
    z = alldat.(season).FixedParams.z;
    zEdges = alldat.(season).FixedParams.zw;
    zWidth = alldat.(season).FixedParams.zwidth; 
    
    Phy_index = alldat.(season).FixedParams.phytoplankton;
    C_index = strcmp(alldat.(season).FixedParams.PP_nut, 'C');
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
        
        
        % all phy C uptake 
        mod = alldat.(season).auxVars.uptake(Phy_index, :, C_index, :, wm_index);
        
        % sum over size class (1st dim)
        mod = squeeze(sum(mod, 1, 'omitnan')); % mmol C /m3 / d
        % integrate over depth until maxdepth
        mod = mod(1:maxdepth, :, :) .* zWidth(1:maxdepth);
        mod = squeeze(sum(mod, 1, 'omitnan'));

        % export for further analyses
        subset = [season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm'];
        productivity.(subset) = mod; 
        
        % convert to g C /m2/d
        mod = mod .* (12.011 * 1e-3);

        % get coordinates
        lat = alldat.(season).Forc.y(:, wm_index);
        lon = alldat.(season).Forc.x(:, wm_index);

        % prepare colorbar ticks
        ll = 0.01;
        %ul = max(mod, [], 'all');
        ul = 3.1;
        Ticks = logspace(log10(ll),log10(ul), 5); 
        TickLabels = round(Ticks, 2, 'significant');
        
        % plot maps:
        

        nexttile
        newmap(proj)
        %g = geobasemap('grayland')
        
        for itraj = 1:size(lat, 2)
            
            trajtable = table(lat(:,itraj), lon(:,itraj), mod(:,itraj), VariableNames=["lat" "lon" "mod"]);
            trajgeotable = table2geotable(trajtable);
            p = geoplot(trajgeotable, ColorVariable='mod', MarkerFaceAlpha=0.5, MarkerSize=0.01);
            if itraj == 1; hold on; end
        end
        p2 = geoplot(land(2:end, :)) % remove antarctica, else the projection produces a bug
        p2.FaceColor = [0.5 0.5 0.5]
        geolimits([65 85], [-20 20])
        hold off
        title([season ' set-up, ' wm ' trajectories'])
        % cb = colorbar('Ticks', Ticks, 'TickLabels', TickLabels);
        % cb.Label.String = 'NPP (g C m^{-2} d^{-1})';
        caxis([ll ul+0.1])
        set(gca,'ColorScale','log');
         
    end

end
colormap(flip(mycolormap))

% styling
t.Padding = 'compact'
t.TileSpacing = 'compact'

% common colorbar
Ticks = [0.01, 0.1, 1, 3]
cb = colorbar('Ticks', Ticks);  % , 'TickLabels', TickLabels
cb.Layout.Tile = "south"
cb.Label.String = 'NPP (g C m^{-2} d^{-1})';


set(gcf, 'Position', [-899  -172   697   464]) % gleiche grosse auch fur export map nehmen
% saveas(t, [Directories.plotDir 'Fig9_NPP_int' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm_allTrajs_waterMass_map.png'])
% hasslich! work-around: aus matlab GUI speichern
clearvars -except autumn summer Directories modTag pSetName productivity




%% export

% make maps of Carbon Export 

maxdepth =8;

watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]

alldat.summer = summer;
alldat.autumn = autumn; 

% map settings
land = readgeotable("landareas.shp");
proj = projcrs(3995, Authority="EPSG"); % 'stereo',  EPSG code https://epsg.io/3995

mycolormap = [165 0 38;
    215 48 39;
    244 109 67;
    253 174 97;
    254 224 144;
    224 243 248;
    171 217 233;
    116 173 209;
    69 117 180;
    49 54 149]./255;



t = tiledlayout(2,2)
title(t, ['Carbon export at depth (' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])

for j = 1:length(seasons)

    season = seasons{j};
    
    C_index = strcmp(alldat.(season).FixedParams.PP_nut, 'C'); % ist eh immer 1
    zwidth = alldat.(season).FixedParams.zwidth(maxdepth); 
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
        
        %Phy
        concPhy = alldat.(season).out.P(:, maxdepth, C_index, :, wm_index);  
        concPhy = squeeze(sum(concPhy, 1, 'omitnan')); % sum over size classes
        %Zoo
        concZoo = alldat.(season).out.Z(:, maxdepth, C_index, :, wm_index);  
        concZoo = squeeze(sum(concZoo, 1, 'omitnan')); % sum over size classes
        %Plankton
        concP = concPhy + concZoo; 
        %OM
        concOM = squeeze(alldat.(season).out.OM(:, maxdepth, C_index, :, wm_index)); 
        
        % extract sinking speeds from FixedParams
        sinksPlankton = alldat.(season).Params.wp(maxdepth, 1); % m/d: careful, if wp is not 0 or different for each size class, then per-sizeclass export has to be calculated before adding 
        sinksOM = [alldat.(season).Params.wDOM1; alldat.(season).Params.wPOM1];

        exOM = concOM .* sinksOM;
        exOM = squeeze(sum(exOM, 1)); % add POM and DOM export (DOM is 0 anyways)

        exP = concP .* sinksPlankton; 

         % get changes in conc due to diffusion
        diffP = alldat.(season).auxVars.B_diffuse(:, maxdepth, C_index, :, wm_index);  % all plankton
        diffP = squeeze(sum(diffP, 1, 'omitnan')) .* zwidth;  % sum over size classes (Phy+Zoo) and muliply with depth layer width

        diffOM = alldat.(season).auxVars.OM_diffuse(:, maxdepth, C_index, :, wm_index); 
        diffOM = squeeze(sum(diffOM, 1, 'omitnan')) .* zwidth; % sum over types (DOM+POM) and multiply with depth layer width

        extot = exOM + exP; % total export of C in mmol m^-2 d^-1 % for all trajs
        % extot = exOM' + exP; % for single trajs
        extot2 = extot + diffOM + diffP; % add changes due to diffusion to changes due to sinking


        subset = [season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm'];
        export.(subset) = extot2; 

        % convert to g C /m2/d
        mod = extot .* (12.011 * 1e-3);
        mod2 = extot2 .* (12.011 * 1e-3);
        
        % get coordinates
        lat = alldat.(season).Forc.y(:, wm_index);
        lon = alldat.(season).Forc.x(:, wm_index);
        
        % plot maps:
        % prepare colorbar ticks
        ll = 0.01;
        %ul = max(mod, [], 'all');
        ul = 0.68;
        Ticks = logspace(log10(ll),log10(ul), 5); 
        TickLabels = round(Ticks, 2, 'significant');



        nexttile
        m = newmap(proj)

         for itraj = 1:size(lat, 2)
            
            trajtable = table(lat(:,itraj), lon(:,itraj), mod2(:,itraj), VariableNames=["lat" "lon" "mod"]);
            trajgeotable = table2geotable(trajtable);
            p = geoplot(trajgeotable, ColorVariable='mod', MarkerFaceAlpha=0.5, MarkerSize=0.01);
            if itraj == 1; hold on; end
        end
        p2 = geoplot(land(2:end, :)) % remove antarctica, else the projection produces a bug
       % p2.FaceColor = [0 0.4470 0.7410]
        p2.FaceColor = [0.5 0.5 0.5]
        geolimits([65 85], [-20 20])
        hold off
        title([season ' set-up, ' wm ' trajectories'])
        % cb = colorbar('Ticks', Ticks, 'TickLabels', TickLabels);
        % cb.Label.String = 'carbon export (g C m^{-2} d^{-1})';
        caxis([ll ul])
        set(gca,'ColorScale','log');

    end

end
colormap(flip(mycolormap))
set(gcf, 'Position', [-899  -172   697   464])

% styling
t.Padding = 'compact'
t.TileSpacing = 'compact'

% common colorbar
Ticks = [0.01, 0.03, 0.1, 0.3, 0.68]
cb = colorbar('Ticks', Ticks);  % , 'TickLabels', TickLabels  
cb.Layout.Tile = "south"
cb.Label.String = 'carbon export (g C m^{-2} d^{-1})';

% and save 
% saveas(t, [Directories.plotDir 'Fig10_export_depthLayer' num2str(maxdepth) '_allTrajs_waterMass_map.png'])

clearvars -except autumn summer Directories modTag pSetName export productivity

% %% Fig 9: plot production and export
% % shows nicely the relationship between the two 
% % as time series
% 
% watermasses = [{'Arctic'}, {'Atlantic'}]
% seasons = [{'summer'}, {'autumn'}]
% 
% 
% alldat.summer = summer;
% alldat.autumn = autumn; 
% 
% 
% maxdepth = 8;
% cmap = lines(2);
% 
% t = tiledlayout(2,2)
% t.Padding = 'compact'
% t.TileSpacing = 'compact'
% %title(t, ['Carbon production and export at depth (' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])
% 
% for j = 1:length(seasons)
% 
%     season = seasons{j};
%     data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
% 
% 
%     for i = 1:length(watermasses)
%         wm = watermasses{i}; 
% 
%         % pick data
%         ex = export.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
%         prod = productivity.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
% 
% 
%         if strcmp(wm, "Arctic"); line_opacity = 0.5; else; line_opacity = 0.1; end
% 
%         % plot
%         nexttile
%         p1 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), prod, 'Color', [cmap(1,:) line_opacity])  %  for multiple trajs
%         hold on
%         p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ex, 'Color', [cmap(2,:) line_opacity])
%         p4 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(prod, 2), 'Color', cmap(1,:), 'LineWidth', 1.5)
%         p5 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(ex, 2), 'Color', cmap(2,:), 'LineWidth', 1.5)
%         ylim([0 250])
%         patch([data_dates fliplr(data_dates)], [min(ylim)*[1 1] max(ylim)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
%             'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
%         hold off
% 
%         %xlabel('time')
%         ylabel('(mmol C m^{-2} d^{-1})')
%         title([season ' setup, ' wm ' trajectories'])
%         legend([p4 p5], {'production', 'export'})
%         set(gca,'children',flipud(get(gca,'children')))
% 
%     end
% end
% linkaxes
% 
% %leg = legend([p4 p5], {'production', 'export'})
% %leg.Position.Tile = 'south'; 
% 
% set(gcf, 'Position', [4   372   808   605])
% 
% saveas(t, [Directories.plotDir 'Fig9_ProdEx' num2str(maxdepth) '_meanallTrajs_waterMass.png'])
% 
% clearvars -except autumn summer Directories modTag pSetName productivity export
% 


%% Fig 9: plot production and export WITH SECONDARY X AXIS
% shows nicely the relationship between the two 
% as time series

watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]


alldat.summer = summer;
alldat.autumn = autumn; 


maxdepth = 8;
cmap = lines(2);

fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 16])
fontsize(fig, 10, 'points')

% t = tiledlayout(2,2)
% t.Padding = 'compact'
% t.TileSpacing = 'compact'
%title(t, ['Carbon production and export at depth (' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])
sp=0;

for j = 1:length(seasons)

    season = seasons{j};
    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];


    referenceDates = data_dates;
    dates = datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum');
    dates_l = datetime("2018-01-01"):datetime("2018-12-31")
    dateDiff = NaN(size(dates_l));
    for k=1:length(dateDiff)
        dat = dates_l(k);
        if dat <= min(referenceDates)
            datDist = days(dat - min(referenceDates));
        elseif dat >= max(referenceDates)
            datDist = days(dat - max(referenceDates));
        else
            datDist = 0;
        end
        dateDiff(k) = datDist;
    end
    dat_i = find(day(dates_l) == 1);
 
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        ex = export.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        prod = productivity.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        
        
        if strcmp(wm, "Arctic"); line_opacity = 0.5; else; line_opacity = 0.1; end

        % plot
        %nexttile
        sp = sp+1;
        subplot(2,2,sp)
        Pos{sp} = get(gca, "Position");
        
        patch([data_dates fliplr(data_dates)], [0*[1 1] 250*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.8,'edgealpha',1)
        
        hold on
        p1 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), prod, 'Color', [cmap(1,:) line_opacity])  %  for multiple traj
        p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), ex, 'Color', [cmap(2,:) line_opacity])
       
        plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(prod, 2), 'Color', [1 1 1], 'LineWidth', 3) % white outline
        p4 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(prod, 2), 'Color', cmap(1,:), 'LineWidth', 1.5)
        plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(ex, 2), 'Color', [1 1 1], 'LineWidth', 3) % white outline
        p5 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(ex, 2), 'Color', cmap(2,:), 'LineWidth', 1.5)
        
        ylim([0 250])
        
        hold off
        
        %xlabel('time')
        ylabel('(mmol C m^{-2} d^{-1})')
       % title([season ' setup, ' wm ' trajectories'])
        % legend([p4 p5], {'production', 'export'})
        % set(gca,'children',flipud(get(gca,'children')))

        xtickformat("MMMMM")
        % secondary x axis
        ax1 = gca;
        ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'YTick', []);
    
        xline(ax2, dates_l(dat_i), 'color', 'none');
        %xline(ax2, dates_l(dat_i), 'color', 'k');
        linkaxes([ax1, ax2]);

         % Customize the secondary x-axis labels
        xticks(ax2, dates_l(dat_i));
        xticklabels(ax2, dateDiff(dat_i))
        xlabel(ax2, 'Temporal distance from target region (d)')
        ax2.XAxis.FontSize = 8;
        % ax2.XLabel.FontSize = 9;
        % add minor labels on the 15th of each month, without label.
        set(ax1,'XMinorTick','on')
        xAx1 = get(ax1,'XAxis');
        xAx1.MinorTickValues=dates_l(dat_i)+15;

        text(10, 225, {[season ' setup,'], [wm ' trajectories']}, 'FontSize',10, 'FontWeight','bold')


    end
end
linkaxes


leg = legend([p4 p5], {'production', 'export'})
leg.Position = [0.405 0.01 0.1683 0.0248] % set position manually
leg.Orientation = "horizontal"


%leg = legend([p4 p5], {'production', 'export'})
%leg.Position.Tile = 'south'; 

% set(gcf, 'Position', [4   372   808   605])

% saveas(gcf, [Directories.plotDir 'Fig9_ProdEx' num2str(maxdepth) '_meanallTrajs_waterMass.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig9_ProdEx' num2str(maxdepth) '_meanallTrajs_waterMass_rev1.png']
saveas(fig, savepath)

clearvars -except autumn summer Directories modTag pSetName productivity export

%% Fig 10 POC maps with binned NPP and export next to it
% create a new figure based on anjas idea
% map with POC or Chl along trajectories. 
% for all ensembles.
% next to it: number for productivity and export fluxes, mean and sd, in 2.5 degree latitude intervals






watermasses = [{'Arctic'}, {'Atlantic'}];
seasons = [{'summer'}, {'autumn'}];

alldat.summer = summer;
alldat.autumn = autumn; 

% map settings
land = readgeotable("landareas.shp");
proj = projcrs(3995, Authority="EPSG"); % 'stereo',  EPSG code https://epsg.io/3995

% ticks for colorbar
TickLabels = [10, 20, 50, 100, 200, 500, 1000]; %0.5, 1, 2, 5, 
Ticks = log10(TickLabels);

% set maxdepth for integration depth
maxdepth = 8; % 100 m 

cmap = colormap("turbo");

fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 18, 18])
fontsize(fig, 10, 'points')

t = tiledlayout(2,2);
%title(t, ['Model trajectories'])
maxval = []

for j = 1:length(seasons)

    season = seasons{j};
    zWidth = alldat.(season).FixedParams.zwidth;
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 

        % pick data
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
             
        
        % get coordinates
        lat = alldat.(season).Forc.y(:, wm_index);
        lon = alldat.(season).Forc.x(:, wm_index);

        % get POC _________________________________________________________
        % det om
        detPOM = squeeze(alldat.(season).out.OM(alldat.(season).FixedParams.POM_index, :, alldat.(season).FixedParams.OM_C_index, :, wm_index));
        
        % Phy POC
        phyvarindex = strcmp(alldat.(season).FixedParams.PP_nut, 'C');
        phyPOM = alldat.(season).out.P(:, :, phyvarindex, :, wm_index);
        % sum over size classes
        phyPOM = squeeze(sum(phyPOM, 1));
       
        
        % Zoo POM
        zoovarindex = strcmp(alldat.(season).FixedParams.ZP_nut, 'C');
        zooPOM = alldat.(season).out.Z(:, :, zoovarindex, :, wm_index);
        % sum over size classes
        zooPOM = squeeze(sum(zooPOM, 1));
       
        
        POC = detPOM + phyPOM + zooPOM; % mmol C /m3

        % integrate POC over depth until maxdepth
        POC_int = POC(1:maxdepth, :, :) .* zWidth(1:maxdepth);
        POC_int = squeeze(sum(POC_int, 1, 'omitnan')); % mmol C /m2

        maxval = [maxval, max(POC_int, [], 'all')];

        % get integrated production of C __________________________________
        phy_index = alldat.(season).FixedParams.phytoplankton;
        prod_C = alldat.(season).auxVars.uptake(phy_index, :, phyvarindex, :, wm_index);
        % sum over size class (1st dim)
        prod_C = squeeze(sum(prod_C, 1, 'omitnan')); % mmol C /m3 / d
        % integrate over depth until maxdepth
        prod_C = prod_C(1:maxdepth, :, :) .* zWidth(1:maxdepth);
        prod_C = squeeze(sum(prod_C, 1, 'omitnan')); % mmol C /m2 / d


        % get export at depth for C _______________________________________
        concP = squeeze(phyPOM(maxdepth, :, :) + zooPOM(maxdepth, :, :)) ; % mmol C /m3
        concOM = squeeze(alldat.(season).out.OM(:, maxdepth, alldat.(season).FixedParams.OM_C_index, :, wm_index)); % mmol C /m3
        
        % extract sinking speeds from FixedParams
        sinksPlankton = alldat.(season).Params.wp(maxdepth, 1); % m/d: careful, if wp is not 0 or different for each size class, then per-sizeclass export has to be calculated before adding 
        sinksOM = [alldat.(season).Params.wDOM1; alldat.(season).Params.wPOM1];
    
        exOM = concOM .* sinksOM;
        exOM = squeeze(sum(exOM, 1)); % add POM and DOM export (DOM is 0 anyways)
    
        exP = concP .* sinksPlankton;

        % get changes in conc due to diffusion
        diffP = alldat.(season).auxVars.B_diffuse(:, maxdepth, phyvarindex, :, wm_index);  % all plankton
        diffP = squeeze(sum(diffP, 1, 'omitnan')) .* zWidth(maxdepth);  % sum over size classes (Phy+Zoo) and muliply with depth layer width
    
        diffOM = alldat.(season).auxVars.OM_diffuse(:, maxdepth, phyvarindex, :, wm_index); 
        diffOM = squeeze(sum(diffOM, 1, 'omitnan')) .* zWidth(maxdepth); % sum over types (DOM+POM) and multiply with depth layer width
    
        extot = exOM + exP; % total export of C in mmol m^-2 d^-1 % for all trajs
        extot2 = extot + diffOM + diffP; % add changes due to diffusion to changes due to sinking

        
        
        % get sampling stations
        scalar_stations = unique(table(alldat.(season).Data.scalar.Latitude, ...
            alldat.(season).Data.scalar.Longitude, VariableNames=["lat", "lon"]));



        % for the side labels: calculate POC, production and export within
        % 2.5 deg latitude bins. (mean of all trajs within each bin)

        minLat = 65; % minimum lat in map and calculation bins
        maxLat = 82.5;  %max lat
        binLat = 2.5;   % step size for map and calculation

        bin_lower = minLat:binLat:maxLat-binLat;
        bin_upper = bin_lower + binLat;
        bin_table = table(bin_lower',bin_upper', NaN(7,1), NaN(7,1), NaN(7,1), 'VariableNames', {'bin_lower', 'bin_upper', 'POC', 'prod', 'ex'});
        % extract data within each bin
        
        for b=1:(maxLat-minLat)/binLat
        
            bin_l = bin_table.bin_lower(b); 
            bin_u = bin_table.bin_upper(b);

            bin_filter = lat > bin_l & lat < bin_u; % filter latitude positions (of each day) that are within this bin
            

            if any(bin_filter, 'all')
                % extract POC within this bin
                bin_table.POC(b) = mean(POC_int(bin_filter)); % mmol C /m2
                bin_table.prod(b) = mean(prod_C(bin_filter)); % mmol C /m2 / d
                bin_table.ex(b) = mean(extot2(bin_filter)); % mmol C /m2 / d

                bin_table.POC_sd(b) = std(POC_int(bin_filter)); % mmol C /m2
                bin_table.prod_sd(b) = std(prod_C(bin_filter)); % mmol C /m2 / d
                bin_table.ex_sd(b) = std(extot2(bin_filter)); % mmol C /m2 / d
             
            end

        end





        
        % plot maps:
      % 
      %   nexttile
      %   m = newmap(proj)
      % 
      %    for itraj = 1:size(lat, 2)
      % 
      %       trajtable = table(lat(:,itraj), lon(:,itraj), log10(POC_int(:,itraj)), VariableNames=["lat" "lon" "POC"]);
      %       trajgeotable = table2geotable(trajtable);
      %       p = geoplot(trajgeotable, ColorVariable='POC', MarkerFaceAlpha=0.5, MarkerSize=0.01);
      %       if itraj == 1; hold on; end
      %    end
      %    % add land
      %   p2 = geoplot(land(2:end, :)) % remove antarctica, else the projection produces a bug
      %  % p2.FaceColor = [0 0.4470 0.7410]
      %   p2.FaceColor = [0.5 0.5 0.5]
      % % add scalar sampling stations
      %   stationgeotable = table2geotable(scalar_stations); 
      %   pStations = geoplot(stationgeotable, "+", MarkerEdgeAlpha=0.3, MarkerSize=3, MarkerEdgeColor="k")
      %   geolimits([65 85], [-20 20])
      %   hold off
      %   title([season ' setup, ' wm ' trajectories'])
      %   %cb = colorbar('Ticks', Ticks, 'TickLabels', TickLabels);
        
        nexttile
        % new try with eqaconic projection
        
        %figure
        worldmap([minLat maxLat],[-20 20])
        ttl = title([season ' setup, ' wm ' trajectories'])
        tpos = ttl.Position
        ttl.Position = [tpos(1) tpos(2)+1.5e5 tpos(3)] % reposition title, so it doesn't overlap with the lat labels
        ax = gca;
        % getm(ax, "MapProjection")
        % show land
        geoshow(land, "FaceColor", [0.5 0.5 0.5])
        hold on
        % show POC along trajectories
        for itraj = 1:size(lat, 2)
            scatterm(lat(:,itraj), lon(:,itraj), 10, log10(POC_int(:,itraj)),'Marker', ".", 'MarkerFaceAlpha', 0.6)
        end
        % show stations
        scatterm(scalar_stations.lat, scalar_stations.lon, 20, 'k', '+')

        % style the maps
        setm(ax, 'PLabelRound', -1);
        setm(ax, 'MLabelLocation', 10);  % Longitude label location
        setm(ax, 'PLabelLocation', binLat);  % Latitude label location
        setm(ax, 'PLineLocation', binLat);  % Grid line location
        
        % setm(ax, 'glinestyle', '-'); % solid gridlines
        setm(ax, 'glinewidth', 0.75);  %    thicker gridlines
        setm(ax, 'flinewidth', 1)  % thinner outter box
        
        
        
        % set limits for color
        clim([log10(10) log10(1600)])

        % Add text labels on the right side
    %textm(80, 20, [{'POC: XX'},{'production: XX'}, {'export: XX'}], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        for b=1:height(bin_table)
            if ~isnan(bin_table.POC(b))
                % textm(bin_table.bin_lower(b), 20, ... % mit SD
                %     [{sprintf('\\bf%.1f±%.1f\\rm', bin_table.prod(b), bin_table.prod_sd(b))}, ...
                %     {sprintf('%.1f±%.1f', bin_table.ex(b), bin_table.ex_sd(b))}],...
                %     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 11)
                textm(bin_table.bin_lower(b), 21, ...
                    [{sprintf('\\bf%.1f\\rm', bin_table.prod(b))}, ...
                    {sprintf('%.1f', bin_table.ex(b))}],...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 11)
            
            end
        end
        
    end

end
% styling
% set(gcf, 'Position', [4  373   626   640])
colormap("turbo")

% common colorbar
cb = colorbar('Ticks', Ticks, 'TickLabels', TickLabels);
%cb.Layout.Tile = "south"; % this overlaps with the labels of the x axis of
%the lower axesm plots 
cb.Title.String="POC [mmol m^{-2}]"
cb.Location = "southoutside"
cb.Position = ([0.13,0.0680,0.775,0.0188]) % manually set the position that would be assigned by cb.Layout.Tile = "south";
cb.Position = cb.Position - [0, 0.04, 0, 0]; % slightly move the colorbar from its original position


 % saveas(t, [Directories.plotDir 'Fig10_POCProdEx' num2str(maxdepth) '_allTrajs.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig10_POCProdEx' num2str(maxdepth) '_allTrajs_rev1.png']
saveas(fig, savepath)

clearvars -except autumn summer Directories modTag pSetName productivity export


%% Fig 11
%% factors shaping changes in size composition
% influence of light limitation, nutrient limitation, and grazing

% 4 panel plot: time vs I_lim, gammaN, V_C:V_N,  and combined grazing losses
% for the example scenario summer/arctic 

vars = {'I_lim', 'gammaN', 'grazing'};
phyIndex = summer.FixedParams.phytoplankton; 

waterMasses = {'Arctic', 'Atlantic'}; 

maxdepth = 6; 
zlayers = summer.FixedParams.z(1:maxdepth);
zwidth = summer.FixedParams.zwidth(1:maxdepth);

sizeClasses = summer.FixedParams.PPdia; 
sizeClassEdges = summer.FixedParams.PPdia_intervals; 

phycolors = {'#cde0b0'; 
    '#95a97b';
    '#5e6b39';
    '#969637';
    '#cdda71';
    '#52933a';
    '#73d670';
    '#d8db33';
    '#7bdb3e'}

% convert colors from hex to rgb
for i=1:length(phycolors)
    str = phycolors{i};
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    phycolors{i} = color;
end


cmap = jet(9);
clab = cell(length(cmap),1)
for s=1:length(sizeClasses)
    clab{s} = [num2str(round(sizeClassEdges(s),2)) ' - ' num2str(round(sizeClassEdges(s+1),2))]; % ' µm'
end


fig = figure
set(fig, 'Units', 'centimeters')
pos = get(fig, 'Position')
set(fig, 'Position', [0, 0, 12, 14])
fontsize(fig, 10, 'points')

t = tiledlayout(length(vars), length(waterMasses),"TileSpacing","compact", "Padding", "compact")
%title(t, {'control of phytoplankton', '(summer set-up)'})



    
for v=1:length(vars)
    
    for wm=1:length(waterMasses)

        waterMass = waterMasses{wm};
        watermassIndex = strcmp(summer.Forc.waterMass, waterMass); 

        var = vars{v};
        nexttile
        switch var 
            case 'I_lim'
                % influence of light limitation
                ilim = summer.auxVars.I_lim(phyIndex, 1:maxdepth, :, watermassIndex); 
                
                % average over trajs
                ilim = mean(ilim, 4);
                % average over depth (weighed)
                ilim = squeeze(sum((ilim .* zwidth' / sum(zwidth)),2));
                
                
                hold on
                for s=1:length(sizeClasses)
                    %clab = [num2str(round(sizeClassEdges(s),2)) ' - ' ...
                    % num2str(round(sizeClassEdges(s+1),2)) ' µm'];
                    plot(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum'), ...
                        ilim(s, :), 'Color', cmap(s,:)) % , 'DisplayName', clab
                end
                hold off
                title({['light limitation'],[waterMass]})
                ylim([0 1])
                ylabel('\gamma L')
            case 'gammaN'
                % gammaN, nutrient limitation
                % it is not saved in the output bus is computed as:
                % Nutrient limitation
                
                N_index = summer.FixedParams.PP_N_index;
                C_index = summer.FixedParams.PP_C_index;
    
                QN = summer.out.P(:,:,N_index,:,:) ./ summer.out.P(:,:,C_index,:,:);
                Qmin_QC = summer.Params.Qmin_QC(phyIndex);
                delQ_QC = summer.Params.delQ_QC(phyIndex);
                
                gammaN = max(0, min(1, (QN - Qmin_QC) ./ delQ_QC));
                
                % pick depth and trajs
                gammaN = gammaN(:,1:maxdepth, :, :, watermassIndex); 
                % average over trajs
                gammaN = mean(gammaN, 5);
                % average over depth
                gammaN = squeeze(sum(gammaN .* zwidth' / sum(zwidth),2))
                
                hold on
                for s=1:length(sizeClasses)
                    %clab = [num2str(round(sizeClassEdges(s),2)) ' - ' num2str(round(sizeClassEdges(s+1),2)) ' µm'];
                    plot(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum'), gammaN(s, :), 'Color', cmap(s,:)) %, 'DisplayName', clab
                end
                hold off
                title({['nutrient limitation'], [waterMass]})
                ylim([0.7 1.01])
                ylabel('\gamma N')
            case 'V_C:V_N'
                % pick data        
                % V_C for phy
                VC = summer.auxVars.V(phyIndex, :, C_index, :, watermassIndex);
                
                % average over depth
                VC = VC(:,1:maxdepth,:,:,:) .* zwidth(1:maxdepth)';
                VC = squeeze(sum(VC./sum(zwidth(1:maxdepth)), 2, 'omitnan'));
                % size(VC)
                VC = mean(VC, 3); % average over trajs
                
                % V_N for phy
                VN = summer.auxVars.V(phyIndex, :, N_index, :, watermassIndex);
                
                % average over depth
                VN = VN(:,1:maxdepth,:,:,:) .* zwidth(1:maxdepth)';
                VN = squeeze(sum(VN./sum(zwidth(1:maxdepth)), 2, 'omitnan'));
                % size(VC)
                VN = mean(VN, 3)
                
                % get C:N ratio
                VCN = VC ./ VN;
                
                hold on
                for s=1:length(sizeClasses)
                    plot(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum'), VCN(s, :), 'Color', cmap(s,:)) 
                end
                hold off
                title('V_C:V_N')
                ylim([0 30])
                ylabel('V_C:V_N')
    
            case 'grazing'
                % data in auxVars.predation_losses
                
                size(summer.auxVars.predation_losses)
                % [ZooP AllP depth nut doy traj]
                C_index = summer.FixedParams.PP_C_index;
                
                predLoss = summer.auxVars.predation_losses(:, phyIndex, 1:maxdepth, C_index, :, watermassIndex); 
                size(predLoss)
                
                % sum over predators -> how much is consumed from each size class
                predLoss = sum(predLoss, 1);
                
                % average over trajs
                predLoss = squeeze(mean(predLoss, 6)); 
                % average over depth
                predLoss = squeeze(sum(predLoss .* zwidth' / sum(zwidth),2)); 
                
                area(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum'), predLoss')
                colororder(cmap)

                
                title({['grazing pressure on phy'], [waterMass]})
                ylabel({'combined grazing','[mmol C m^{-3} d^{-1}]'})

                

        end % var switch
    xtickformat("MMMMM");
    ax = gca;
    ax.YMinorTick = "on"
    end % end var
end % end wm

leg = legend(clab, 'NumColumns', 3)
title(leg,'size classes [µm]') 
leg.Layout.Tile = 'south';

% set(gcf, 'Position', [-944  -260   448   639])
% saveas(t, [Directories.plotDir 'Fig11_phy_lim' num2str(round(sum(summer.FixedParams.zwidth(1:maxdepth)))) 'summer.png'])
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/revised_figs/', 'Fig11_phy_lim' num2str(round(sum(summer.FixedParams.zwidth(1:maxdepth)))) 'summer_rev1.png']
saveas(fig, savepath)

%% from here on, left over plots that are either moved to appendix or deleted entirely follow. 
% I have not yet cleaned this up. 


%% pe ratio & f-ratio together

% as time series
watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]


alldat.summer = summer;
alldat.autumn = autumn; 


maxdepth = 8;
cmap = lines(2);
% prodThreshold = 1e-4; % threshold, above which f ratio can be calculated (otherwise produciton is too low)


t = tiledlayout(2,2)
title(t, ['Carbon export efficiency and f-ratio at depth (' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])

for j = 1:length(seasons)

    season = seasons{j};
    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];

    z = alldat.(season).FixedParams.z;
    zEdges = alldat.(season).FixedParams.zw;
    zWidth = alldat.(season).FixedParams.zwidth; 
    
    Phy_index = alldat.(season).FixedParams.phytoplankton; % only Phy takes up DIN
    N_P_index = strcmp(alldat.(season).FixedParams.PP_nut, 'N');
    N_OM_index = strcmp(alldat.(season).FixedParams.OM_nut, 'N');
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);

        % pick pe data
        ex = export.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        prod = productivity.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        
        % calculate movinf average (one sided backwards) for prod and ex to
        % smooth it a bit (loop over trajs)
        ex_smooth = NaN(size(ex)); 
        prod_smooth = NaN(size(prod)); 
        for j=1:size(ex, 2)
            ex_smooth(:,j) = moving_average_backw(ex(:,j), 5);
            prod_smooth(:,j) = moving_average_backw(prod(:,j), 5); 
        end

        ex = ex_smooth; 
        prod = prod_smooth;

        % threshold, above which production happens and calulation of pe
        % makes sense
        threshold = prod < 0.6625*10; % evtl hoch oder runter setzen 
        
        ex2 = ex;
        ex2(threshold) = 0;
      

        % ex3 = round(ex, 3)
        % prod3 = round(prod,3)

        pe = ex2 ./ prod;
       % pe(threshold) = NaN;
        
        % f-ratio
        % all phy N uptake 
        VNtot = alldat.(season).auxVars.uptake(Phy_index, :, N_P_index, :, wm_index); % N uptake mmol N /m3 / d
        
        % sum over size class (1st dim)
        VNtot = squeeze(sum(VNtot, 1, 'omitnan')); 
        % integrate over depth until maxdepth
        VNtot = VNtot(1:maxdepth, :, :) .* zWidth(1:maxdepth);
        VNtot = squeeze(sum(VNtot, 1, 'omitnan')); % mmol N /m2 / d

        % all N remineralisation
        VNremin = alldat.(season).auxVars.OM_remin(:, :, N_OM_index, :, wm_index); 
        % sum over OM types (add DON and PON)
        VNremin = squeeze(sum(VNremin, 1, 'omitnan')); % mmol N /m3 / d
         % integrate over depth until maxdepth
        VNremin = VNremin(1:maxdepth, :,:) .* zWidth(1:maxdepth);
        VNremin = squeeze(sum(VNremin, 1, 'omitnan')); % mmol N /m2 / d


        % calculate moving average (one sided backwards) for VNtot and VNremin to
        % smooth it a bit (loop over trajs)
        VNremin_smooth = NaN(size(VNremin)); 
        VNtot_smooth = NaN(size(VNtot)); 
        for j=1:size(ex, 2)
            VNremin_smooth(:,j) = moving_average_backw(VNremin(:,j), 5);
            VNtot_smooth(:,j) = moving_average_backw(VNtot(:,j), 5); 
        end

        VNremin = VNremin_smooth;
        VNtot = VNtot_smooth; 

        % % get N productivity again, as a threshold for f-ratio calc. 
        % prodN = alldat.(season).auxVars.uptake(Phy_index, :, N_P_index, :, wm_index);
        % % sum over size class (1st dim)
        % prodN = squeeze(sum(prodN, 1, 'omitnan')); % mmol N /m3 / d
        % % integrate over depth until maxdepth
        % prodN = prodN(1:maxdepth, :, :) .* zWidth(1:maxdepth);
        % prodN = squeeze(sum(prodN, 1, 'omitnan')); % mmol N /m2 / d
        
        % calculate f-ratio

        % set VNremin = 0, wenn VNtotal < VNremin
        filter = VNtot < VNremin + 1; % 0.01 mmol /m3 /d weil uber 100m int.
        VNremin_zero = VNremin;
        VNremin_zero(filter) = 0
        fratio = (VNtot - VNremin_zero) ./ VNtot; 

        

      %  fratio(VNtot <= 0.01) = NaN; 

      % time average prior to actual date berechnen 
      % prod ex VNtot VNremin (nicht erst pe und f-ratio)
      % dann vllt thresholds wieder verringern


        % plot
        nexttile
        p1 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), pe, 'Color', [cmap(1,:) 0.1], 'LineWidth', 0.3)
        hold on
        p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), fratio, 'Color', [cmap(2,:) 0.1], 'LineWidth', 0.3)
        p3 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(pe,2), 'Color', [cmap(1,:) 1], 'LineWidth', 1)
        p4 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(fratio,2), 'Color', [cmap(2,:) 1], 'LineWidth', 1)
        xlabel('time')
        ylabel('pe-ratio / f-ratio')
        ylim([-0.1 1.1])
        patch([data_dates fliplr(data_dates)], [min(ylim)*[1 1] max(ylim)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.2,'edgealpha',0.5)
        title([season ' set-up, ' wm ' trajectories'])
        legend([p3(1) p4(1)], {'pe', 'f-ratio'})
        hold off

    end
end
linkaxes
ylim([-0.1 1.1])

t.TileSpacing = "compact"
t.Padding = "compact"

set(gcf, 'Position', [4   372   808   605])

saveas(t, [Directories.plotDir 'Fig12_pe_efficiency_f_ratio' num2str(maxdepth) '_allTrajs_waterMass.png'])

clearvars -except autumn summer Directories modTag pSetName productivity export

%% plot pe for C and N, instead of f-ratio

% load data generated in plot_NPP_export_pe_efficiency.m
productivity_N = load('productivity_N.mat').productivity_N
export_N = load('export_N.mat').export_N

% as time series
watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]


alldat.summer = summer;
alldat.autumn = autumn; 


maxdepth = 8;
cmap = lines(2);
% prodThreshold = 1e-4; % threshold, above which f ratio can be calculated (otherwise produciton is too low)


t = tiledlayout(2,2)
title(t, ['Export efficiency of Carbon and Nitrogen at depth (' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) ' m)'])

for j = 1:length(seasons)

    season = seasons{j};
    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];

    z = alldat.(season).FixedParams.z;
    zEdges = alldat.(season).FixedParams.zw;
    zWidth = alldat.(season).FixedParams.zwidth; 
    
    Phy_index = alldat.(season).FixedParams.phytoplankton; % only Phy takes up DIN
    N_P_index = strcmp(alldat.(season).FixedParams.PP_nut, 'N');
    N_OM_index = strcmp(alldat.(season).FixedParams.OM_nut, 'N');
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 
        wm_index = strcmp(alldat.(season).Forc.waterMass, wm);

        % pick pe data
        ex_C = export.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        prod_C = productivity.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        ex_N = export_N.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        prod_N = productivity_N.([season '_' wm '_' num2str(abs(alldat.summer.FixedParams.zw(maxdepth+1))) 'm']);
        
        % % calculate movinf average (one sided backwards) for prod and ex to
        % % smooth it a bit (loop over trajs)
        % ex_smooth = NaN(size(ex)); 
        % prod_smooth = NaN(size(prod)); 
        % for j=1:size(ex, 2)
        %     ex_smooth(:,j) = moving_average_backw(ex(:,j), 5);
        %     prod_smooth(:,j) = moving_average_backw(prod(:,j), 5); 
        % end
        % 
        % ex = ex_smooth; 
        % prod = prod_smooth;

        % threshold, above which production happens and calulation of pe
        % makes sense
        threshold = prod_C < 0.6625*10; % evtl hoch oder runter setzen 
        %threshold = prod_C < 0.6625*10 & ex_C < 0.6625*10;
        % 
        % ex2 = ex;
        % ex2(threshold) = 0;
        % 
        % 
        % % ex3 = round(ex, 3)
        % % prod3 = round(prod,3)

        % ex_C(threshold) = 0;
        % ex_N(threshold) = 0;

        pe_C = ex_C ./ prod_C;
        pe_N = ex_N ./ prod_N;
       
        pe_C(threshold) = NaN;
        pe_N(threshold) = NaN;
       

        % plot
        nexttile
        p1 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), pe_C, 'Color', [cmap(1,:) 0.1], 'LineWidth', 0.3)
        hold on
        p2 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), pe_N, 'Color', [cmap(2,:) 0.1], 'LineWidth', 0.3)
        p3 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(pe_C,2, "omitnan"), 'Color', [cmap(1,:) 1], 'LineWidth', 1)
        p4 = plot(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'), mean(pe_N,2, "omitnan"), 'Color', [cmap(2,:) 1], 'LineWidth', 1)
        xlabel('time')
        ylabel('pe-ratio')
        ylim([0 1.1])
        patch([data_dates fliplr(data_dates)], [min(ylim)*[1 1] max(ylim)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.2,'edgealpha',0.5)
        title([season ' set-up, ' wm ' trajectories'])
        legend([p3(1) p4(1)], {'Carbon', 'Nitrogen'})
        hold off
        

    end
end
linkaxes
ylim([0 1.1])

t.TileSpacing = "compact"
t.Padding = "compact"

set(gcf, 'Position', [4   372   808   605])

saveas(t, [Directories.plotDir 'Fig12_pe_C_N_' num2str(maxdepth) '_allTrajs_waterMass.png'])

clearvars -except autumn summer Directories modTag pSetName productivity_N export_N productivity export


%% variability pdfs, depth averaged

% DIN, Chl, POC and DOC
% during sampling time

% plot diffKDE of state vars, for sampling period, discriminated arctic and
% atlantic trajectories by line type/color
% depth averaged



alldat.summer = summer;
alldat.autumn = autumn;

vars = [{'DIN'}, {'Chl'}, {'POC'}]; % DIN bekommt eigene Kachel, {'DOC'}
seasons = [{'summer'}, {'autumn'}];
waterMasses = [{'Arctic'}, {'Atlantic'}]; 
maxdepth = 6; % up to 62m

figure
t = tiledlayout(length(seasons), length(vars));
title(t, ['density estimates of tracer concentrations (mean over upper ' num2str(round(sum(summer.FixedParams.zwidth(1:maxdepth)))) ' m)' ]);

cmap = lines(length(waterMasses));

for s = 1:length(seasons)
    season = seasons{s};

    modDates = datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum');
    modMonths = month(modDates); 
    zwidth = alldat.(season).FixedParams.zwidth(1:maxdepth);
    
    % pick model equivalents with a wider time frame
    samplingDates = unique(alldat.(season).Data.scalar.Yearday); 
    
        
    % var loop for individual tiles
    for v = 1:length(vars)
        var = vars{v};
        nexttile
        hold on
    
        % wm loop for linetypes
        for wm = 1:length(waterMasses)
            waterMass = waterMasses{wm};
            wmIndex = strcmp(alldat.(season).Forc.waterMass, waterMass);


            % pick data
            switch var
                case 'DIN'
                    mod = squeeze(alldat.(season).out.N(:,1:maxdepth,samplingDates, wmIndex)); % mmol N/m3
                    xlimits = [0 16];
                    xlab = 'mmol N m^{-3}';
                    xmax = max(alldat.(season).out.N(:,1:maxdepth,:, :), [], 'all'); 
                    xmax = ceil(xmax*1.1);
                case 'Chl'
                    ChlIndex = alldat.(season).FixedParams.PP_Chl_index;
                    mod = alldat.(season).out.P(:,1:maxdepth, ChlIndex, samplingDates, wmIndex);  % mg Chl/m3
                    % sum over size classes (1st dimension)
                    mod = squeeze(sum(mod, 1));
                    xlab = 'mg Chl m^{-3}';
                    xlimits = [0 3]
                    xmax = max(sum(alldat.(season).out.P(:,1:maxdepth, ChlIndex, :, :),1), [], 'all'); 
                    xmax = ceil(xmax*1.1);
                case 'POC'
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, samplingDates, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_C_index = alldat.(season).FixedParams.PP_C_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_C_index, samplingDates, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_C_index = alldat.(season).FixedParams.ZP_C_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, samplingDates, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all POC
                    mod = modDet + modPhy + modZoo;
                    xlab = 'mmol C m^{-3}';
                    xlimits = [0 35];

                    xmax = max(squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, :, :)) ...
                        + squeeze(sum(alldat.(season).out.P(:,1:maxdepth, PP_C_index, :, :),1)) ...
                        + squeeze(sum(alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, :, :),1)), ...
                        [], 'all')
                    xmax = ceil(xmax)
                case 'DOC'
                    % detritus
                    DOMindex = alldat.(season).FixedParams.DOM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    mod = squeeze(alldat.(season).out.OM(DOMindex,1:maxdepth, OM_C_index, samplingDates, wmIndex));  % mmol C/m3
                    
                    xlab = 'mmol C m^{-3}';
                    xlimits = [0 15];

                    xmax = max(squeeze(alldat.(season).out.OM(DOMindex,1:maxdepth, OM_C_index, :, :)),[], 'all')
                    xmax = ceil(xmax);

            end

            % average over depth (gewichtet nach Tiefe)
            mod = squeeze(sum(mod .* zwidth, 1)./ sum(zwidth)); 

            % reshape to long vector
            mod = reshape(mod, [1, numel(mod)]);

            xmin = 0; 
            
            % calculate KDE
            try
                KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax); % sometimes diffKDE fails because of numerical problems in the discretisation
            catch
                warning('oh oh')
                % try
                %     KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax+0.5*xmax);
                %     warning([char(datetime(1,m,1, 'Format', 'MMMM')) ', ' waterMass ' trajectories, ' var ': used exaggerated xmin and xmax due to numerical problems'])
                % catch
                %     KDE = py.diffKDE.KDE(py.numpy.array(mod));
                %     warning([char(datetime(1,m,1, 'Format', 'MMMM')) ', ' waterMass ' trajectories, ' var ': used default xmin and xmax due to numerical problems'])
                % end
            end

            KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
            KDE_y = double(KDE{1}); % convert first numpy array to matlab double
            KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double

             % fix the ugly hanging ends in the air..
             KDE_x_t = [KDE_x(1), KDE_x, KDE_x(end)];
             KDE_y_t = [1e-6, KDE_y, 1e-6];

             plot(KDE_x, KDE_y, 'Color', cmap(wm,:))

           %  if wm == 1; hold on;  end; %else hold off;

     %   end

        % und beobachtungen hinzufugen

        switch var
            case 'DIN'
                dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'N');
            case 'Chl'
                dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'chl_a');
            case 'POC'
                dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'POC');
                   
        end

        % create a watermass vector for filtering Data (from events and watermasses) 
        wm_obs = alldat.(season).Data.scalar.waterMass(alldat.(season).Data.scalar.Event); 
        % alldat.(season).Data.scalar.Event is 865x1
        % alldat.(season).Data.scalar.Event is 36x1

        % get indeces for watermass of observations.  
        % in summer, do it 'exact', but in autumn there are no 'pure'
        % observations, therefore consider the mixed ones as "arctic" 
        if(strcmp(season, "summer"))
            if strcmp(waterMass, "Arctic"); wm_obs_index = strcmp(wm_obs, "Arctic"); else; wm_obs_index = strcmp(wm_obs, "Atlantic"); end
        else
            if strcmp(waterMass, "Arctic"); wm_obs_index = strcmp(wm_obs, "Arctic")|strcmp(wm_obs, "Arctic/Atlantic"); else; wm_obs_index = strcmp(wm_obs, "Atlantic"); end
        end
        
        % get data for current var
        dat = alldat.(season).Data.scalar.Value(dat_var_index&wm_obs_index); % mmol /m3 or mg /m3

        % only continue, if there is data for this var in this setup.
        if isempty(dat)

            
            leg = {['model - ' waterMasses{1}],  ['model - ' waterMasses{2}]};
            % plot legend 
        else

            % get depth of data
            data_depth = alldat.(season).Data.scalar.Depth(dat_var_index&wm_obs_index);
            % get events of data
            data_events = alldat.(season).Data.scalar.Event(dat_var_index&wm_obs_index);
            data_events_u = unique(data_events); 
    
            % interpolate to 1m steps
            data_interp = [];
            depth_interp = [];
            for i_event = 1:length(data_events_u)
                % extract this event, as all events have different depths
                ev = data_events_u(i_event);
                dat_ev = dat(data_events == ev);
                dep_ev = data_depth(data_events == ev);
    
                if length(dat_ev) == 1
                    dat_ev_interp = dat_ev;
                else
    
    
                    % interpolate
                    
                    % create target depth vector
                    dep_ev_interp = 1:dep_ev(end);
    
                    % generate very small numbers to add to dep_ev, incase some
                    % depths are non-unique (this leads to failing of
                    % interpolation
                    rn = rand([length(dep_ev) 1])/100;
    
                    % interp1
                    dat_ev_interp = interp1(dep_ev+rn, dat_ev,dep_ev_interp);
        
                    % figure
                    % plot(dat_ev, -dep_ev, "o")
                    % hold on
                    % plot(dat_ev_interp, -dep_ev_interp)
    
                    % retain only shallower than 62m
                    dat_ev_interp = dat_ev_interp(dep_ev_interp < 63);
                    dep_ev_interp = dep_ev_interp(dep_ev_interp < 63);
        
                    % get water column average
                    dat_ev_interp = mean(dat_ev_interp, 'omitnan');
                end
       
                % append to output vectors
                data_interp = [data_interp, dat_ev_interp];
                %depth_interp = [depth_interp, dep_ev_interp];
    
            end
    
    
            % exclude nans, if not done already during averaging
            data_interp = data_interp(~isnan(data_interp));
            
            xmax = ceil(max(data_interp))*1.2;
            % calc diffKDE 
            KDE = py.diffKDE.KDE(py.numpy.array(data_interp), xmin=0, xmax=xmax); % sometimes diffKDE fails because of numerical problems in the discretisation
            
            KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
            KDE_y = double(KDE{1}); % convert first numpy array to matlab double
            KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double
    
             % fix the ugly hanging ends in the air..
             KDE_x_t = [KDE_x(1), KDE_x, KDE_x(end)];
             KDE_y_t = [1e-6, KDE_y, 1e-6];
            
             % plot
             plot(KDE_x, KDE_y, 'LineStyle',':', 'Color', cmap(wm,:))
            
            

            leg = {['model - ' waterMasses{1}], ['observations - ' waterMasses{1}], ['model - ' waterMasses{2}], ['observations - ' waterMasses{2}]}
        end

        end
        hold off
         % hier plot schoen machen
        legend(leg, 'Location', 'best')
        xlabel(xlab)
        ylabel('prob. density')
        title([season, ', ' var]);
        xlim(xlimits)
       % ylim([0 1])
    
        
    end
    
end

t.TileSpacing = "compact"
t.Padding = "compact"

set(gcf, 'Position', [-1318        -136         926         515])
saveas(t, [Directories.plotDir 'Fig13a_diffKDE_PDF_conc_depth_avg_separate_stateVars_samplingTime_with_obs.png'])

clearvars -except autumn summer Directories modTag pSetName
%% variability pdfs

% DIN, Chl, POC and DOC
% during sampling time

% plot diffKDE of state vars, for sampling period, discriminated arctic and
% atlantic trajectories by line type/color
% depth resolved



alldat.summer = summer;
alldat.autumn = autumn;

vars = [{'DIN'}, {'Chl'}, {'POC'}, {'DOC'}]; % DIN bekommt eigene Kachel
seasons = [{'summer'}, {'autumn'}];
waterMasses = [{'Arctic'}, {'Atlantic'}]; 
maxdepth = 6; % up to 62m

figure
t = tiledlayout(length(seasons), length(vars));
title(t, ['density estimates of tracer concentrations (resolved over upper ' num2str(round(sum(summer.FixedParams.zwidth(1:maxdepth)))) ' m)' ]);

cmap = lines(length(waterMasses));

for s = 1:length(seasons)
    season = seasons{s};

    modDates = datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum');
    modMonths = month(modDates); 
    zwidth = alldat.(season).FixedParams.zwidth(1:maxdepth);
    
    % pick model equivalents with a wider time frame
    samplingDates = unique(alldat.(season).Data.scalar.Yearday); 
    
        
    % var loop for individual tiles
    for v = 1:length(vars)
        var = vars{v};
        nexttile
    
        % wm loop for linetypes
        for wm = 1:length(waterMasses)
            waterMass = waterMasses{wm};
            wmIndex = strcmp(alldat.(season).Forc.waterMass, waterMass);


            % pick data
            switch var
                case 'DIN'
                    mod = squeeze(alldat.(season).out.N(:,1:maxdepth,samplingDates, wmIndex)); % mmol N/m3
                    xlimits = [0 15];
                    xlab = 'mmol N m^{-3}';
                    xmax = max(alldat.(season).out.N(:,1:maxdepth,:, :), [], 'all'); 
                    xmax = ceil(xmax*1.1);
                case 'Chl'
                    ChlIndex = alldat.(season).FixedParams.PP_Chl_index;
                    mod = alldat.(season).out.P(:,1:maxdepth, ChlIndex, samplingDates, wmIndex);  % mg Chl/m3
                    % sum over size classes (1st dimension)
                    mod = squeeze(sum(mod, 1));
                    xlab = 'mg Chl m^{-3}';
                    xlimits = [0 4]
                    xmax = max(sum(alldat.(season).out.P(:,1:maxdepth, ChlIndex, :, :),1), [], 'all'); 
                    xmax = ceil(xmax*1.1);
                case 'POC'
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, samplingDates, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_C_index = alldat.(season).FixedParams.PP_C_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_C_index, samplingDates, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_C_index = alldat.(season).FixedParams.ZP_C_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, samplingDates, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all POC
                    mod = modDet + modPhy + modZoo;
                    xlab = 'mmol C m^{-3}';
                    xlimits = [0 22];

                    xmax = max(squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, :, :)) ...
                        + squeeze(sum(alldat.(season).out.P(:,1:maxdepth, PP_C_index, :, :),1)) ...
                        + squeeze(sum(alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, :, :),1)), ...
                        [], 'all')
                    xmax = ceil(xmax)
                case 'DOC'
                    % detritus
                    DOMindex = alldat.(season).FixedParams.DOM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    mod = squeeze(alldat.(season).out.OM(DOMindex,1:maxdepth, OM_C_index, samplingDates, wmIndex));  % mmol C/m3
                    
                    xlab = 'mmol C m^{-3}';
                    xlimits = [0 30];

                    xmax = max(squeeze(alldat.(season).out.OM(DOMindex,1:maxdepth, OM_C_index, :, :)),[], 'all')
                    xmax = ceil(xmax);

            end

            % average over depth (gewichtet nach Tiefe)
            %mod = squeeze(sum(mod .* zwidth, 1)./ sum(zwidth)); 
            % OR: resolve over depth
            depthsMOD = abs(alldat.(season).FixedParams.z(1:maxdepth));
            depthsTARGET = 1:sum(zwidth);
  
            modINT = interpolate_mod(mod, depthsMOD, depthsTARGET'); % interpolate over depth
            size(modINT)

            % reassign modINT to mod
            mod = modINT; 

            % reshape to long vector
            mod = reshape(mod, [1, numel(mod)]);

            xmin = 0; 
            
            % calculate KDE
            try
                KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax); % sometimes diffKDE fails because of numerical problems in the discretisation
            catch
                warning('oh oh')
                % try
                %     KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax+0.5*xmax);
                %     warning([char(datetime(1,m,1, 'Format', 'MMMM')) ', ' waterMass ' trajectories, ' var ': used exaggerated xmin and xmax due to numerical problems'])
                % catch
                %     KDE = py.diffKDE.KDE(py.numpy.array(mod));
                %     warning([char(datetime(1,m,1, 'Format', 'MMMM')) ', ' waterMass ' trajectories, ' var ': used default xmin and xmax due to numerical problems'])
                % end
            end

            KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
            KDE_y = double(KDE{1}); % convert first numpy array to matlab double
            KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double

             % fix the ugly hanging ends in the air..
             KDE_x_t = [KDE_x(1), KDE_x, KDE_x(end)];
             KDE_y_t = [1e-6, KDE_y, 1e-6];

             plot(KDE_x, KDE_y, 'Color', cmap(wm,:))

             if wm == 1; hold on; else hold off; end;

        end
        % hier plot schoen machen
        leg = legend(waterMasses)
        xlabel(xlab)
        ylabel('prob. density')
        title([season, ', ' var]);
        xlim(xlimits)
       % ylim([0 1])
    
    end
    
end

t.TileSpacing = "compact"
t.Padding = "compact" 

set(gcf, 'Position', [-1318        -136         926         515])
saveas(t, [Directories.plotDir 'Fig13b_diffKDE_PDF_conc_depth_resolved_separate_stateVars_samplingTime.png'])

clearvars -except autumn summer Directories modTag pSetName

%% DOC pdfs over time
% by months (not two month groups) 

% DOC over time, nur sommer
% distinguish between water masses, each month
% use diffKDE instead
% averagoe over upper water column

depthLevel_index = 1:6; % upper 62 m
depthLevelEdges = [summer.FixedParams.zw(1) summer.FixedParams.zw(5+1)] ;
watermasses = [{'Arctic'}, {'Atlantic'}]


% https://www.learnui.design/tools/data-color-picker.html#divergent
colorPalette = ["#003f5c", "#444e86", "#955196", "#dd5182", "#ff6e54", "#ffa600"]
%colorPalette = colormap(parula(6));

figure
t = tiledlayout(1,2)
title(t, 'Density estimates of DOC concentrations (upper 62m), summer set-up') %average over 


% get indeces for OM
DOM_index = summer.FixedParams.DOM_index;
OM_C_index = summer.FixedParams.OM_C_index;

% convert time (datenum) to months
months = month(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum')); 
u_months = unique(months);
monthsNames = month(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum'), 'shortname'); 
u_monthsNames = unique(monthsNames, 'stable');  

for wm = 1:length(watermasses)
    watermass = watermasses{wm}; 

    wm_index = strcmp(summer.Forc.waterMass, watermass);
    zwidth = summer.FixedParams.zwidth(depthLevel_index); 
    % begin plotting
    nexttile
    hold on
    for i = 5:length(u_months) % loop through months, start in may
        
        % set time index
        time_index = months == i;

        % pick model output
        mod = squeeze(summer.out.OM(DOM_index, depthLevel_index, OM_C_index, time_index, wm_index));
        size(mod)

        % average over depth
        % mod = sum(mod .* zwidth, 1)./ sum(zwidth)

        % OR: resolve over depth
        depthsMOD = abs(summer.FixedParams.z(depthLevel_index));
        depthsTARGET = 1:sum(zwidth);
  
        modINT = interpolate_mod(mod, depthsMOD, depthsTARGET'); % interpolate over depth
        size(modINT)

        % reassign modINT to mod
        mod = modINT; 

        % reshape to long vector
        mod = reshape(mod, [1, numel(mod)]);
       
        % get KDE
        KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=0, xmax=30) % diffKDE
        KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
        KDE_y = double(KDE{1}); % convert first numpy array to matlab double
        KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double

        % fix the ugly hanging ends in the air..
        KDE_x_t = [KDE_x(1), KDE_x, KDE_x(end)];
        KDE_y_t = [1e-6, KDE_y, 1e-6];
        
        % add KDE to plot
        if(i == 2); hold on; end
        depthDisp = u_monthsNames{i};
        plot(KDE_x_t, KDE_y_t, 'DisplayName', depthDisp, 'Color', colorPalette(i-4)) % diffKDE
        
    end
    hold off
    title([watermass ' trajectories'])
    ylabel('probability density')
    xlabel('[DOC] (mmol m^{-3})')
     xlim([0 30])
     ylim([0 0.4]) 
end
leg = legend('NumColumns',6)
leg.Layout.Tile = 'south';
linkaxes

t.TileSpacing = "compact"
t.Padding = "compact"

set(gcf, 'Position', [680   614   520   264])
saveas(t, [Directories.plotDir 'Fig15_diffKDE_DOC_over_time_summer_watermass.png'])

clearvars -except autumn summer Directories modTag pSetName



%% calculate trajectory distances

dist_summer = abs((summer.Forc.y(1,:)-summer.Forc.y(end,:)))

[m, I] = min(dist_summer)
m
max(dist_summer)
size(summer.Forc.y(1,:))

plot(summer.Forc.x(:,I), summer.Forc.y(:,I))
summer.Forc.waterMass(I) % the small min is because this trajectory recirculates and comes back to same lat where it started

arc = strcmp(summer.Forc.waterMass, "Arctic")
min(abs((summer.Forc.y(1,arc)-summer.Forc.y(end,arc))))
max(abs((summer.Forc.y(1,arc)-summer.Forc.y(end,arc))))

arc = strcmp(autumn.Forc.waterMass, "Arctic")
min(abs((autumn.Forc.y(1,arc)-autumn.Forc.y(end,arc))))
max(abs((autumn.Forc.y(1,arc)-autumn.Forc.y(end,arc))))


dist_autumn = abs((autumn.Forc.y(1,:)-autumn.Forc.y(end,:)))

[m, I] = min(dist_autumn)
m
max(dist_autumn)