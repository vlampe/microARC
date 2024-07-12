% create a table with statisitics for the state variables during
% observation periods
% DIN, Chl, POC, DOC; averaged over the upper 100 m
% summer arc, summer atl, autumn arc, autumn atl 

close all
clear all
clearvars -except autumn summer Directories modTag pSetName



modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_fullRunForPaper';
%% add paths 
addpath(genpath('~/Documents/microARC model')) % location of my work with the model
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

%% 

str = struct()
str.mean = 1
str.mode = 2
str.median = 3

tab = struct2table(str)
tab2 = tab.*2
[tab; tab2]

m.key = 'id' % var, season, wm
m.mode = 3
m.median = 3
m.mean = 3
m.CI_l = 2
m.CI_u =4

m = struct2table(m)

 tab = array2table(zeros(0,6), 'VariableNames',{'key', 'mode', 'median', 'mean', 'CI_l', 'CI_u'})

[tab; m;m]

%%

resTab = array2table(zeros(0,7), 'VariableNames',{'key', 'unit', 'mode', 'median', 'mean', 'CI_l', 'CI_u'})

alldat.summer = summer;
alldat.autumn = autumn; 

maxdepth = 8; % 100m, end of euphotic zone (from visual analysis of PhyC/PhyChl profiles)

watermasses = [{'Arctic'}, {'Atlantic'}];
seasons = [{'summer'}, {'autumn'}];
vars = [{'DIN'}, {'Chl'}, {'POC'}, {'PON'}, {'POC:PON'}, {'DOC'}]; 

timeFrame = 0; % include simulations +- 7 days around sampling of observations

for v=1:length(vars)
    var = vars{v};
    var_row = cell2table({string(var) NaN NaN NaN NaN NaN NaN}, 'VariableNames', ...
        {'key', 'unit', 'mode', 'median', 'mean', 'CI_l', 'CI_u'}); 
     resTab = [resTab ;var_row];
    for s=1:length(seasons)
        season = seasons{s};
    
        % pick model equivalents with a wider time frame
        samplingDates = unique(alldat.(season).Data.scalar.Yearday); 
        timeIndex = samplingDates(1)-timeFrame:samplingDates(end)+timeFrame; 

        zwidth = alldat.(season).FixedParams.zwidth(1:maxdepth);
        % implement calculation of variability params for observations
        %
            % pick data according to var
            switch var
                case 'DIN'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'N');
                    datUnit = 'mmol N m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index); % mmol /m3 or mg Chl /m3
                case 'Chl'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'chl_a');
                    datUnit = 'mg Chl a m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index); % mmol /m3 or mg Chl /m3
                case 'POC'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'POC');
                    datUnit = 'mmol C m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index); % mmol /m3 or mg Chl /m3
                case 'PON'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'PON');
                    datUnit = 'mmol N m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index); % mmol /m3 or mg Chl /m3
                case 'POC:PON'
                    datUnit = 'mol:mol';
                    dat_C_index = strcmp(alldat.(season).Data.scalar.Variable, 'POC');
                    dat_var_index = dat_C_index;
                    datC = alldat.(season).Data.scalar.Value(dat_C_index); % mmol /m3 or mg Chl /m3
                    dat_N_index = strcmp(alldat.(season).Data.scalar.Variable, 'PON');
                    datN = alldat.(season).Data.scalar.Value(dat_N_index);  
                    
                    if isempty(datN)
                        dat = [];
                    else
                        dat = datC ./ datN; % mol:mol
                    end
                case 'DOC'
                    dat = [];
                    datUnit = 'mmol C m^{-3}';
            end
            

             if isempty(dat)

                % create empty labelled row
                %m.key = string([var ' ' season ' observations']); 
                m.key = string([season ' observations']); 
                m.unit = {datUnit};
                m.mode = NaN; % get mode from KDE instead. When each number is only found once (likely with the high decimal precision), mode function returns the first one...
                m.median = NaN; 
                m.mean = NaN;
                m.CI_l = NaN;
                m.CI_u = NaN;
            else
                % interpolate, average and get KDE

                % get depth of data
                data_depth = alldat.(season).Data.scalar.Depth(dat_var_index);
                % get events of data
                data_events = alldat.(season).Data.scalar.Event(dat_var_index);
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
                        dat_ev_interp = dat_ev_interp(dep_ev_interp <= 100);
                        dep_ev_interp = dep_ev_interp(dep_ev_interp <= 100);
            
                        % get water column average
                        dat_ev_interp = mean(dat_ev_interp, 'omitnan');
                    end
           
                    % append to output vectors
                    data_interp = [data_interp, dat_ev_interp];
                    %depth_interp = [depth_interp, dep_ev_interp];
        
                end
           
            % get KDE and statistics
            
            % exclude nans, if not done already during averaging
            data_interp = data_interp(~isnan(data_interp));
            
            xmax = ceil(max(data_interp))*1.2;
            % calc diffKDE 
            KDE = py.diffKDE.KDE(py.numpy.array(data_interp), xmin=0, xmax=xmax); % sometimes diffKDE fails because of numerical problems in the discretisation
            
            KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
            KDE_y = double(KDE{1}); % convert first numpy array to matlab double
            KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double
    
            
             % % plot
             % plot(KDE_x, KDE_y, 'Color', 'k')
            
            h = diff(KDE_x(1:2)); 
            pdf = KDE_y;

            cdf = cumsum(pdf)*h; 
           
            % figure
            % plot(KDE_x, pdf)
            % hold on
            % plot(KDE_x, cdf)

            % find confidence intervals
            % in y vector of cdf, find position where y is closest to .025,
            % (0.5 for median), and .975
            [temp, pos] = min(abs(cdf - 0.025));
            m.CI_l = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.5));
            m.median = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.975));
            m.CI_u = KDE_x(pos);

            % find mode from KDE: where is KDE max?
            [temp, pos] = max(pdf);
            m.mode = KDE_x(pos);
            
            m.mean = mean(data_interp); 
            %m.key = string([var ' ' season ' observations']); 
            m.key = string([season ' observations']); 
            m.unit = {datUnit};

            

             end
            % append row to table
            % create row for table
            % convert m to table 
            m = struct2table(m);   
            resTab = [resTab; m];
            %vertcat(resTab, m)'

            clear m KDE* pdf cdf  


%
    
        for w=1:length(watermasses)
            watermass = watermasses{w};
            wmIndex = strcmp(alldat.(season).Forc.waterMass, watermass);
    
            % disp([var '_' season '_' watermass])
            % ok reihenfolge stimmt schon mal
            
            % pick model output
            switch var
                case 'DIN'
                    mod = squeeze(alldat.(season).out.N(:,1:maxdepth,timeIndex, wmIndex)); % mmol N/m3
                    % convert to mg C / m3
                    %mod = mod .* (12.011);
                    %modUnit = "mg C m^{-3}";
                    modUnit = "mmol C m^{-3}";
                case 'Chl'
                    ChlIndex = alldat.(season).FixedParams.PP_Chl_index;
                    mod = alldat.(season).out.P(:,1:maxdepth, ChlIndex, timeIndex, wmIndex);  % mg Chl/m3
                    % sum over size classes (1st dimension)
                    mod = squeeze(sum(mod, 1));
                    modUnit = "mg Chl a m^{-3}";
                case 'POC'
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_C_index = alldat.(season).FixedParams.PP_C_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_C_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_C_index = alldat.(season).FixedParams.ZP_C_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all POC
                    mod = modDet + modPhy + modZoo;
                    % convert to mg C / m3
                    %mod = mod .* (12.011);
                    %modUnit = "mg C m^{-3}";
                    modUnit = "mmol C m^{-3}";

                case 'PON'
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_N_index = alldat.(season).FixedParams.OM_N_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_N_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_N_index = alldat.(season).FixedParams.PP_N_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_N_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_N_index = alldat.(season).FixedParams.ZP_N_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_N_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all PON
                    mod = modDet + modPhy + modZoo;
                    % convert to mg C / m3
                    %mod = mod .* (14.007);
                    %modUnit = "mg N m^{-3}";
                    modUnit = "mmol N m^{-3}";
                case 'POC:PON'
                    %  CARBON
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_C_index = alldat.(season).FixedParams.PP_C_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_C_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_C_index = alldat.(season).FixedParams.ZP_C_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all POC
                    modC = modDet + modPhy + modZoo;

                    % NITROGEN
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_N_index = alldat.(season).FixedParams.OM_N_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_N_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_N_index = alldat.(season).FixedParams.PP_N_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_N_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_N_index = alldat.(season).FixedParams.ZP_N_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_N_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all PON
                    modN = modDet + modPhy + modZoo;
                    
                    % get C:N
                    mod = modC./modN;
                    modUnit = "POC:PON [mol:mol]";

                case 'DOC'
                    DOMindex = alldat.(season).FixedParams.DOM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    mod = squeeze(alldat.(season).out.OM(DOMindex,1:maxdepth, OM_C_index, timeIndex, wmIndex)); 
                    % convert to mg C / m3
                    % mod = mod .* (12.011);
                    % modUnit = "mg C m^{-3}";
                    modUnit = "mmol C m^{-3}";
            end

            % average over depth (gewichtet nach Tiefe)
            mod = squeeze(sum(mod .* zwidth, 1)./ sum(zwidth));

            % reshape to long vector
            mod = reshape(mod, [1, numel(mod)]);

            % statistiken berechnen
            m.key = string([season ' model: ' watermass]); 
            %m.key = string([var ' ' season ' model: ' watermass]); 
            m.unit = modUnit;
            % m.mode = mode(mod);
            m.mode = NaN; % get mode from KDE instead. When each number is only found once (likely with the high decimal precision), mode function returns the first one...
            % m.median = median(mod);
            m.median = NaN; 
            m.mean = mean(mod);
            
            
            m.CI_l = NaN;
            m.CI_u = NaN;

            xmax = ceil(max(mod) + 0.1 * range(mod)); % oder 0.1, wenn KDE nicht in der Luft hangen darf
            xmin = 0;
            % xmax = ceil(max(mod));
            % xmin = 0; 


            % KDE pdf berechnen, dann CDF, dann CIs rausfinden
            % KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax);
             try
                    KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax); % sometimes diffKDE fails because of numerical problems in the discretisation
             catch
                    %warning('no KDE computed due to numerical problems')
                    try
                        KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin-0.5*xmin, xmax=xmax+0.5*xmax);
                        warning([char(m.key) ': used exaggerated xmin and xmax due to numerical problems'])
                    catch
                        KDE = py.diffKDE.KDE(py.numpy.array(mod));
                        warning([char(m.key) ': used default xmin and xmax due to numerical problems'])
                    end
            end

            KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
            KDE_y = double(KDE{1}); % convert first numpy array to matlab double
            KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double
            h = diff(KDE_x(1:2)); 
            pdf = KDE_y;

            cdf = cumsum(pdf)*h; 
           
            % figure
            % plot(KDE_x, pdf)
            % hold on
            % plot(KDE_x, cdf)

            % find confidence intervals
            % in y vector of cdf, find position where y is closest to .025,
            % (0.5 for median), and .975
            [temp, pos] = min(abs(cdf - 0.025));
            m.CI_l = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.5));
            m.median = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.975));
            m.CI_u = KDE_x(pos);

            % find mode from KDE: where is KDE max?
            [temp, pos] = max(pdf);
            m.mode = KDE_x(pos);
            
            % convert m to table 
            m = struct2table(m); 

            % append to result table
            resTab = [resTab; m];
            %vertcat(resTab, m)'

            clear m KDE* pdf cdf  
    
        end
    end
end
resTab2 = resTab;
resTab2.dCI = resTab2.CI_u - resTab2.CI_l;

resTab2{:, 3:end} = round(resTab2{:, 3:end},1)



%writetable(resTab2, [Directories.plotDir 'tracerconcentrations.csv'])
%clearvars -except autumn summer Directories modTag pSetName resTab

%% same as above, but here, observations are split by watermasse, too
% although I am not conviinced this is good, as the 'mixed' stations in
% autumn are 90% atlantic and only 10% arctic.....

resTab = array2table(zeros(0,7), 'VariableNames',{'key', 'unit', 'mode', 'median', 'mean', 'CI_l', 'CI_u'})

alldat.summer = summer;
alldat.autumn = autumn; 

maxdepth = 8; % 100m, end of euphotic zone (from visual analysis of PhyC/PhyChl profiles)

watermasses = [{'Arctic'}, {'Atlantic'}];
seasons = [{'summer'}, {'autumn'}];
vars = [{'DIN'}, {'Chl'}, {'POC'}, {'PON'}, {'POC:PON'}, {'DOC'}]; 

timeFrame = 0; % include simulations +- 7 days around sampling of observations

for v=1:length(vars)
    var = vars{v};
    var_row = cell2table({string(var) NaN NaN NaN NaN NaN NaN}, 'VariableNames', ...
        {'key', 'unit', 'mode', 'median', 'mean', 'CI_l', 'CI_u'}); 
     resTab = [resTab ;var_row];
    for s=1:length(seasons)
        season = seasons{s};
    
        % pick model equivalents with a wider time frame
        samplingDates = unique(alldat.(season).Data.scalar.Yearday); 
        timeIndex = samplingDates(1)-timeFrame:samplingDates(end)+timeFrame; 

        zwidth = alldat.(season).FixedParams.zwidth(1:maxdepth);



%
    
        for w=1:length(watermasses)
            watermass = watermasses{w};
            wmIndex = strcmp(alldat.(season).Forc.waterMass, watermass);
    
            % pick observations
            % implement calculation of variability params for observations
        
            % create a watermass vector for filtering Data (from events and watermasses) 
            wm_obs = alldat.(season).Data.scalar.waterMass(alldat.(season).Data.scalar.Event); 
            % alldat.(season).Data.scalar.Event is 865x1
            % alldat.(season).Data.scalar.Event is 36x1

            % get indeces for watermass of observations.  
            % in summer, do it 'exact', but in autumn there are no 'pure'
            % observations, therefore consider the mixed ones as "arctic" 
            if(strcmp(season, "summer"))
                if strcmp(watermass, "Arctic"); wm_obs_index = strcmp(wm_obs, "Arctic"); else; wm_obs_index = strcmp(wm_obs, "Atlantic"); end
            else
                if strcmp(watermass, "Arctic"); wm_obs_index = strcmp(wm_obs, "Arctic")|strcmp(wm_obs, "Arctic/Atlantic"); else; wm_obs_index = strcmp(wm_obs, "Atlantic"); end
            end
            

            % pick data according to var
            switch var
                case 'DIN'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'N');
                    datUnit = 'mmol N m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index&wm_obs_index); % mmol /m3 or mg Chl /m3
                case 'Chl'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'chl_a');
                    datUnit = 'mg Chl a m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index&wm_obs_index); % mmol /m3 or mg Chl /m3
                case 'POC'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'POC');
                    datUnit = 'mmol C m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index&wm_obs_index); % mmol /m3 or mg Chl /m3
                case 'PON'
                    dat_var_index = strcmp(alldat.(season).Data.scalar.Variable, 'PON');
                    datUnit = 'mmol N m^{-3}';
                    % get data for current var
                    dat = alldat.(season).Data.scalar.Value(dat_var_index&wm_obs_index); % mmol /m3 or mg Chl /m3
                case 'POC:PON'
                    datUnit = 'mol:mol';
                    dat_C_index = strcmp(alldat.(season).Data.scalar.Variable, 'POC');
                    dat_var_index = dat_C_index;
                    datC = alldat.(season).Data.scalar.Value(dat_C_index&wm_obs_index); % mmol /m3 or mg Chl /m3
                    dat_N_index = strcmp(alldat.(season).Data.scalar.Variable, 'PON');
                    datN = alldat.(season).Data.scalar.Value(dat_N_index&wm_obs_index);  
                    
                    if isempty(datN)
                        dat = [];
                    else
                        dat = datC ./ datN; % mol:mol
                    end
                case 'DOC'
                    dat = [];
                    datUnit = 'mmol C m^{-3}';
            end
            

             if isempty(dat)

                % create empty labelled row
                %m.key = string([var ' ' season ' observations']); 
                m.key = string([season ' observations: ' watermass]); 
                m.unit = {datUnit};
                m.mode = NaN; % get mode from KDE instead. When each number is only found once (likely with the high decimal precision), mode function returns the first one...
                m.median = NaN; 
                m.mean = NaN;
                m.CI_l = NaN;
                m.CI_u = NaN;
            else
                % interpolate, average and get KDE

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
                        dat_ev_interp = dat_ev_interp(dep_ev_interp <= 100);
                        dep_ev_interp = dep_ev_interp(dep_ev_interp <= 100);
            
                        % get water column average
                        dat_ev_interp = mean(dat_ev_interp, 'omitnan');
                    end
           
                    % append to output vectors
                    data_interp = [data_interp, dat_ev_interp];
                    %depth_interp = [depth_interp, dep_ev_interp];
        
                end
           
            % get KDE and statistics
            
            % exclude nans, if not done already during averaging
            data_interp = data_interp(~isnan(data_interp));
            
            xmax = ceil(max(data_interp))*1.2;
            % calc diffKDE 
            KDE = py.diffKDE.KDE(py.numpy.array(data_interp), xmin=0, xmax=xmax); % sometimes diffKDE fails because of numerical problems in the discretisation
            
            KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
            KDE_y = double(KDE{1}); % convert first numpy array to matlab double
            KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double
    
            
             % % plot
             % plot(KDE_x, KDE_y, 'Color', 'k')
            
            h = diff(KDE_x(1:2)); 
            pdf = KDE_y;

            cdf = cumsum(pdf)*h; 
           
            % figure
            % plot(KDE_x, pdf)
            % hold on
            % plot(KDE_x, cdf)

            % find confidence intervals
            % in y vector of cdf, find position where y is closest to .025,
            % (0.5 for median), and .975
            [temp, pos] = min(abs(cdf - 0.025));
            m.CI_l = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.5));
            m.median = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.975));
            m.CI_u = KDE_x(pos);

            % find mode from KDE: where is KDE max?
            [temp, pos] = max(pdf);
            m.mode = KDE_x(pos);
            
            m.mean = mean(data_interp); 
            %m.key = string([var ' ' season ' observations']); 
            m.key = string([season ' observations: ' watermass]); 
            m.unit = {datUnit};

            

             end
            % append row to table
            % create row for table
            % convert m to table 
            m = struct2table(m);   
            resTab = [resTab; m];
            %vertcat(resTab, m)'

            clear m KDE* pdf cdf  







            % _____________________________________________________________
            
            % pick MODEL output
            switch var
                case 'DIN'
                    mod = squeeze(alldat.(season).out.N(:,1:maxdepth,timeIndex, wmIndex)); % mmol N/m3
                    % convert to mg C / m3
                    %mod = mod .* (12.011);
                    %modUnit = "mg C m^{-3}";
                    modUnit = "mmol C m^{-3}";
                case 'Chl'
                    ChlIndex = alldat.(season).FixedParams.PP_Chl_index;
                    mod = alldat.(season).out.P(:,1:maxdepth, ChlIndex, timeIndex, wmIndex);  % mg Chl/m3
                    % sum over size classes (1st dimension)
                    mod = squeeze(sum(mod, 1));
                    modUnit = "mg Chl a m^{-3}";
                case 'POC'
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_C_index = alldat.(season).FixedParams.PP_C_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_C_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_C_index = alldat.(season).FixedParams.ZP_C_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all POC
                    mod = modDet + modPhy + modZoo;
                    % convert to mg C / m3
                    %mod = mod .* (12.011);
                    %modUnit = "mg C m^{-3}";
                    modUnit = "mmol C m^{-3}";

                case 'PON'
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_N_index = alldat.(season).FixedParams.OM_N_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_N_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_N_index = alldat.(season).FixedParams.PP_N_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_N_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_N_index = alldat.(season).FixedParams.ZP_N_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_N_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all PON
                    mod = modDet + modPhy + modZoo;
                    % convert to mg C / m3
                    %mod = mod .* (14.007);
                    %modUnit = "mg N m^{-3}";
                    modUnit = "mmol N m^{-3}";
                case 'POC:PON'
                    %  CARBON
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_C_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_C_index = alldat.(season).FixedParams.PP_C_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_C_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_C_index = alldat.(season).FixedParams.ZP_C_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_C_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all POC
                    modC = modDet + modPhy + modZoo;

                    % NITROGEN
                    % detritus
                    POMindex = alldat.(season).FixedParams.POM_index;
                    OM_N_index = alldat.(season).FixedParams.OM_N_index;
                    modDet = squeeze(alldat.(season).out.OM(POMindex,1:maxdepth, OM_N_index, timeIndex, wmIndex));  % mmol C/m3
                    % phytoplankton
                    PP_N_index = alldat.(season).FixedParams.PP_N_index; 
                    modPhy = alldat.(season).out.P(:,1:maxdepth, PP_N_index, timeIndex, wmIndex);  % mmol C/m3
                    % sum over size classes
                    modPhy = squeeze(sum(modPhy, 1));
                    % zooplankton
                    ZP_N_index = alldat.(season).FixedParams.ZP_N_index; 
                    modZoo = alldat.(season).out.Z(:,1:maxdepth, ZP_N_index, timeIndex, wmIndex);  % moml C/m3
                    % sum over size classes
                    modZoo = squeeze(sum(modZoo, 1));

                    % add all PON
                    modN = modDet + modPhy + modZoo;
                    
                    % get C:N
                    mod = modC./modN;
                    modUnit = "POC:PON [mol:mol]";

                case 'DOC'
                    DOMindex = alldat.(season).FixedParams.DOM_index;
                    OM_C_index = alldat.(season).FixedParams.OM_C_index;
                    mod = squeeze(alldat.(season).out.OM(DOMindex,1:maxdepth, OM_C_index, timeIndex, wmIndex)); 
                    % convert to mg C / m3
                    % mod = mod .* (12.011);
                    % modUnit = "mg C m^{-3}";
                    modUnit = "mmol C m^{-3}";
            end

            % average over depth (gewichtet nach Tiefe)
            mod = squeeze(sum(mod .* zwidth, 1)./ sum(zwidth));

            % reshape to long vector
            mod = reshape(mod, [1, numel(mod)]);

            % statistiken berechnen
            m.key = string([season ' model: ' watermass]); 
            %m.key = string([var ' ' season ' model: ' watermass]); 
            m.unit = modUnit;
            % m.mode = mode(mod);
            m.mode = NaN; % get mode from KDE instead. When each number is only found once (likely with the high decimal precision), mode function returns the first one...
            % m.median = median(mod);
            m.median = NaN; 
            m.mean = mean(mod);
            
            
            m.CI_l = NaN;
            m.CI_u = NaN;

            xmax = ceil(max(mod) + 0.1 * range(mod)); % oder 0.1, wenn KDE nicht in der Luft hangen darf
            xmin = 0;
            % xmax = ceil(max(mod));
            % xmin = 0; 


            % KDE pdf berechnen, dann CDF, dann CIs rausfinden
            % KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax);
             try
                    KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin, xmax=xmax); % sometimes diffKDE fails because of numerical problems in the discretisation
             catch
                    %warning('no KDE computed due to numerical problems')
                    try
                        KDE = py.diffKDE.KDE(py.numpy.array(mod), xmin=xmin-0.5*xmin, xmax=xmax+0.5*xmax);
                        warning([char(m.key) ': used exaggerated xmin and xmax due to numerical problems'])
                    catch
                        KDE = py.diffKDE.KDE(py.numpy.array(mod));
                        warning([char(m.key) ': used default xmin and xmax due to numerical problems'])
                    end
            end

            KDE = cell(KDE); % convert py touple to matlab cell (cell contains 2 numpy arrays)
            KDE_y = double(KDE{1}); % convert first numpy array to matlab double
            KDE_x = double(KDE{2}); % convert 2nd numpy array to matlab double
            h = diff(KDE_x(1:2)); 
            pdf = KDE_y;

            cdf = cumsum(pdf)*h; 
           
            % figure
            % plot(KDE_x, pdf)
            % hold on
            % plot(KDE_x, cdf)

            % find confidence intervals
            % in y vector of cdf, find position where y is closest to .025,
            % (0.5 for median), and .975
            [temp, pos] = min(abs(cdf - 0.025));
            m.CI_l = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.5));
            m.median = KDE_x(pos);

            [temp, pos] = min(abs(cdf - 0.975));
            m.CI_u = KDE_x(pos);

            % find mode from KDE: where is KDE max?
            [temp, pos] = max(pdf);
            m.mode = KDE_x(pos);
            
            % convert m to table 
            m = struct2table(m); 

            % append to result table
            resTab = [resTab; m];
            %vertcat(resTab, m)'

            clear m KDE* pdf cdf  
    
        end
    end
end
resTab2 = resTab;
resTab2.dCI = resTab2.CI_u - resTab2.CI_l;

resTab2{:, 3:end} = round(resTab2{:, 3:end},1)

writetable(resTab2, [Directories.plotDir 'tracerconcentrations_obs_and_mod.csv'])
