% investigate growth rates by sze class

% vars = {'I_lim', 'gammaN', 'grazing'};
vars = {'VN', 'VC'};
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
            case 'VN'
                % growth rate
                N_index = summer.FixedParams.PP_N_index;
                VN = summer.auxVars.V(phyIndex, 1:maxdepth, N_index, :, watermassIndex); 
                
                % average over trajs
                VN = mean(VN, 5);
                % average over depth (weighed)
                VN = squeeze(sum((VN .* zwidth' / sum(zwidth)),2));
                
                
                hold on
                for s=1:length(sizeClasses)
                    %clab = [num2str(round(sizeClassEdges(s),2)) ' - ' ...
                    % num2str(round(sizeClassEdges(s+1),2)) ' µm'];
                    plot(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum'), ...
                        VN(s, :), 'Color', cmap(s,:)) % , 'DisplayName', clab
                    end
                hold off
                title({['V_N'],[waterMass]})
               % ylim([0 1])
                ylabel('V_N (d^{-1})')
                
            case 'VC'
                % growth rate
                C_index = summer.FixedParams.PP_C_index;
                VC = summer.auxVars.V(phyIndex, 1:maxdepth, C_index, :, watermassIndex); 
                
                % average over trajs
                VC = mean(VC, 5);
                % average over depth (weighed)
                VC = squeeze(sum((VC .* zwidth' / sum(zwidth)),2));
                
                
                hold on
                for s=1:length(sizeClasses)
                    %clab = [num2str(round(sizeClassEdges(s),2)) ' - ' ...
                    % num2str(round(sizeClassEdges(s+1),2)) ' µm'];
                    plot(datetime(summer.Forc.t(:,1),'ConvertFrom','datenum'), ...
                        VC(s, :), 'Color', cmap(s,:)) % , 'DisplayName', clab
                    end
                hold off
                title({['V_C'],[waterMass]})
              %  ylim([0 1])
                ylabel('V_C (d^{-1})')

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
                
                gammaN = max(0, min(1, (QN - Qmin_QC) ./ delQ_QC)); % nutrient uptake regulation term Qstat, muss auch noch hoch h
                
                
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

savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/R2_VC_VN_time.png']
saveas(fig, savepath)


clearvars -except autumn summer Directories modTag pSetName productivity export
%%

% plot realised growth rate for different times (diff light and nutrient regimes) over cell size

fig = figure
phyIndex = summer.FixedParams.phytoplankton; 

waterMasses = {'Arctic', 'Atlantic'}; 

maxdepth = 6; 
zlayers = summer.FixedParams.z(1:maxdepth);
zwidth = summer.FixedParams.zwidth(1:maxdepth);

sizeClasses = summer.FixedParams.PPdia; 
sizeClassEdges = summer.FixedParams.PPdia_intervals; 


t = tiledlayout(1, length(waterMasses),"TileSpacing","compact", "Padding", "compact")
%title(t, {'control of phytoplankton', '(summer set-up)'})

for wm=1:length(waterMasses)

    waterMass = waterMasses{wm};
    watermassIndex = strcmp(summer.Forc.waterMass, waterMass); 

    nexttile

    % growth rate
    C_index = summer.FixedParams.PP_C_index;
    VC = summer.auxVars.V(phyIndex, 1:maxdepth, C_index, :, watermassIndex); 
    
    % average over trajs
    VC = mean(VC, 5);
    % average over depth (weighed)
    VC = squeeze(sum((VC .* zwidth' / sum(zwidth)),2));
    
            

    % x: cell size, 9 classes, sizeClasses micro m
    % VC, realised carbon production rate d^-1

    timesteps = 90:21:260  % size(VC, 2)
    cmap = turbo(length(timesteps))

    hold on
    for t=1:length(timesteps)
        plot(sizeClasses, VC(:,timesteps(t)), 'Color', cmap(t,:))
    end
    
    xlim([1,200])
    ax = gca
    set(ax, 'XScale', 'log');
    xlabel('cell size / ESD (µm)')
    ylabel('V_C (d^{-1})')

    ax.XAxis.TickValues = [1 2 5 10 20 50 100 200];
    ax.XTickLabelRotation = 0;
   
    title({['V_C'],[waterMass]})

end % end wm


leg = legend(string(datetime(summer.Forc.t(timesteps,1),'ConvertFrom','datenum')))
leg.Layout.Tile = "south"
leg.Orientation = "vertical"
leg.NumColumns = 3
% 
savepath = ['~/Documents/microARC/Manuscripts/Manuscript microARC Lagrangian modelling plankton variability/REVIEWS/R2_VC_size.png']
saveas(fig, savepath)
