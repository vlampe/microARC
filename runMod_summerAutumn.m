% make a full (or optionally partial) model run for a specific selected
% parameter set
% to produce model output for second paper (summer-autumn transition)

% 
% Refresh the workspace
clear; clc; close all; delete(gcp('nocreate'));

% run all trajectories?
runAllTrajectories = true; 

% how should Params be set? (which method)
setParamMethod = 'method3';

%% set tag for this model run for saving of outputs

% % optimisation by Aidan Hunter for the model with corrected diffusion
% modTag = 'optimisedParams_AH';
% pSetName = 'fittedParameters_RMS_Hellinger_ZPratio_Atlantic_correctedDiffusion.mat'
% % uses fittedParameters_RMS_Hellinger_ZPratio_Atlantic_correctedDiffusion.mat
% 


% use Aidans opt params, but larger pmax_b AND smaller Qmin_QC_a AND larger
% theta. use only icefree trajectories for the arctic in summer and adjust
% faulty DIN values form SINMOD simulation for initial conditions


% some calibration runs:
% modTag = 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin_fullRunForPaper'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_5_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_a_3.5_pmax_b_-0.1_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.04_theta_5.6_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_Gmax_a_11_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.04_theta_5.6_Gmax_a_11_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.035_theta_5.6_Gmax_a_11_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_5_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_0.5_icelessArcSummer_adjustedVin' % almost no grazing
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_22_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_16_icelessArcSummer_adjustedVin'
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin_closed'
% modTag = 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin_fullRunForPaper'
% modTag = 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin_fullRunForPaper'
% modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_adjustedVin_fullRunForPaper'; % for comparisions of variability
%modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_rPOM_0.02_icelessArcSummer_adjustedVin'


% final model config used for the submitted manuscript:
modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_fullRunForPaper'; 

% % for discussion: a run for only summer atlantic, with increased wPOM. 
% modTag = 'optimisedParams_AH_wPOM_100_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_AtlanticOnly_forDiscussion'; 
% modTag = 'optimisedParams_AH_wPOM_100_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_forDiscussion'; 
% 
% % for discussion: a run for only summer arctic, with decreased aP. 
% modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_1.915e-07_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'; 
% 
% % for discussion: a run for only summer arctic, with increased aP. 
% modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_7.66e-07_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'; 
% 
% % for discussion: a run for only summer arctic, with even further increased aP. 
% modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_30e-07_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'; 
% 
% % for discussion: a run for only summer arctic, with even further increased aP. 
% modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_3.83e-06_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'; 

% new trial with increased wPOM: 
%% add paths 
% addpath(genpath('../Documents/microARC model/'));
% directory ensemble_pack now contains needed functions for model set up,
% matching, and data and trajectories for autumn

% Include microARC path
addpath(genpath(fileparts(which('run_model'))));


% set directories
Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', []);
display(Directories)

Directories.resultsDir  = ['../Documents/microARC model/2nd paper/mod output/' modTag '/'];
Directories.plotDir  = ['../Documents/microARC model/2nd paper/mod output/' modTag '/plots/'];

if ~exist(Directories.resultsDir)
    mkdir(Directories.resultsDir)
end
if ~exist(Directories.plotDir)
    mkdir(Directories.plotDir)
end
    


%% set up model
% Use default set-up to load and organise data and initialise model parameters.
% Set-up may be modified here by passing some extra arguments as name-value
% pairs (preferable), or directly modified within modelSetUp.m
% summer
[Forc_s, FixedParams_s, Params_s, Data_s] = modelSetUp2(Directories, ...    
    'displayAllOutputs', true, 'seasonConfig', 'summer', 'years', 2018); 
% autumn
[Forc_a, FixedParams_a, Params_a, Data_a] = modelSetUp2(Directories, ...
    'displayAllOutputs', true, 'seasonConfig', 'autumn', 'years', 2018, ...
    'maxDist',100); % maxDist is 100 here, because when using default (25) all trajs would be of atlantic origin


% map the tjaectories and choose one arctic and one atlantic
map_trajectories(Forc_s, FixedParams_s);
% id 2935 for arc, 621 atl
map_trajectories(Forc_a, FixedParams_a);
% id 291 for arc, 220 atl

%% optional: manipulate Forc to only model fewer trajectories instead of full set

if contains(modTag, 'icelessArcSummer')
    % instead using singular trajs as in all other quick simulations, use
    % an iceless trajectory for summer/arctic --- for testing purposes
    arcticTraj = 2008
end

if runAllTrajectories == false
    % summer
    Forc_s_ = Forc_s;
    fieldnames = fields(Forc_s);

    
    if ~exist('arcticTraj', 'var')
        arcticTraj = 3126;
    end
    
    %summerTrajs = [343 3126];
    summerTrajs = [343 arcticTraj];

    summerTrajsIndex = [find(Forc_s.iTraj == summerTrajs(1)) find(Forc_s.iTraj == summerTrajs(2))];
    nTrajOrig = Forc_s.nTraj;

    for i = 1:length(fieldnames)

        if size(Forc_s.(fieldnames{i}),2) == nTrajOrig
            Forc_s_.(fieldnames{i}) = Forc_s.(fieldnames{i})(:,summerTrajsIndex);
        elseif size(Forc_s.(fieldnames{i}),3) == nTrajOrig
            Forc_s_.(fieldnames{i}) = Forc_s.(fieldnames{i})(:,:,summerTrajsIndex);
        end
    end
    Forc_s_.nTraj = length(summerTrajs);

    % create a list of trajIDs per event, for later use in matchModOutput2Data2
    evTraj_s = Data_s.scalar.evTraj;
    evTrajID_s = Forc_s.iTraj(evTraj_s);

    Forc_s = Forc_s_;
    %autumn
    Forc_a_ = Forc_a;
    fieldnames = fields(Forc_a);

    autumnTrajs = [220 294];
    autumnTrajsIndex = [find(Forc_a.iTraj == autumnTrajs(1)) find(Forc_a.iTraj == autumnTrajs(2))];


    nTrajOrig = Forc_a.nTraj;

    for i = 1:length(fieldnames)

        if size(Forc_a.(fieldnames{i}),2) == nTrajOrig
            Forc_a_.(fieldnames{i}) = Forc_a.(fieldnames{i})(:,autumnTrajsIndex);
        elseif size(Forc_a.(fieldnames{i}),3) == nTrajOrig
            Forc_a_.(fieldnames{i}) = Forc_a.(fieldnames{i})(:,:,autumnTrajsIndex);
        end
    end
    Forc_a_.nTraj = length(autumnTrajs);

    % create a list of trajIDs per event, for later use in matchModOutput2Data2
    evTraj_a = Data_a.scalar.evTraj;
    evTrajID_a = Forc_a.iTraj(evTraj_a);

    Forc_a = Forc_a_;

    clear Forc_s_ Forc_a_ fieldnames nTrajOrig

    % map again to check if selection was successful
    map_trajectories(Forc_s, FixedParams_s);
    map_trajectories(Forc_a, FixedParams_a);

end

% if a full set of trajectories should be simulated, but only ice-free ones
% in summer/arctic:
if (runAllTrajectories == true) && contains(modTag, 'icelessArcSummer')
    % summer
    Forc_s_ = Forc_s;
    fieldnames = fields(Forc_s);
    
    atlanticTrajs = Forc_s.iTraj(strcmp(Forc_s.waterMass, 'Atlantic'));    % all of them
    arcticTrajs = [2008 2033 2015 2010 2131 2002 2011 2019 2181 2018];   % hand picked based on light availability/ice cover


    % for discussion: do some extra runs with just arctic or just atlantic
    % trajs.
    if contains(modTag, 'ArcticOnly')
        atlanticTrajs = []
    elseif contains(modTag, 'AtlanticOnly')
        arcticTrajs = []
    end

   
    summerTrajs = [atlanticTrajs arcticTrajs];

    % summerTrajsIndex = [find(Forc_s.iTraj == summerTrajs(1)) find(Forc_s.iTraj == summerTrajs(2))];
    summerTrajsIndex = find(any(Forc_s_.iTraj == summerTrajs'));
    nTrajOrig = Forc_s.nTraj;

    for i = 1:length(fieldnames)

        if size(Forc_s.(fieldnames{i}),2) == nTrajOrig
            Forc_s_.(fieldnames{i}) = Forc_s.(fieldnames{i})(:,summerTrajsIndex);
        elseif size(Forc_s.(fieldnames{i}),3) == nTrajOrig
            Forc_s_.(fieldnames{i}) = Forc_s.(fieldnames{i})(:,:,summerTrajsIndex);
        end
    end
    Forc_s_.nTraj = length(summerTrajs);

    % create a list of trajIDs per event, for later use in matchModOutput2Data2
    evTraj_s = Data_s.scalar.evTraj;
    evTrajID_s = Forc_s.iTraj(evTraj_s);

    Forc_s = Forc_s_;
end



%% SET PARAMS
% Params are equal for summer and autumnConfig. 

switch setParamMethod
    case 'method1' % method 1: pick one parameter set from the ensemble runs
        % they are the ensemble member 183 
        memberID = 22;   % member 183 looks nice

        Params = Params_s; 
        %update params

        % load parameter value file
        pEnsemble = load([Directories.resultsDir 'parameterEnsemble2000.txt']);

        % and  assign parameter names
        Pnames = { 'A', 'h', 'm2', 'aP', 'theta', 'xi', 'aG', 'sigG', 'Lambda', ...
        'lambda_max','Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a',...
        'Qmax_delQ_b', 'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b',...
        'pmax_a', 'pmax_b', 'Gmax_a', 'Gmax_b', 'm_a', 'm_b', 'wp_a', 'wp_b',...
        'beta1', 'beta2', 'beta3', 'wPOM1', 'rDON' , 'rPON', 'rDOC', 'rPOC'};

        ensArgIn = {};
        for j = 1:length(Pnames)
           ensArgIn{1,j} = Pnames{j};
           ensArgIn{2,j} = pEnsemble(memberID,j);  % memberID with index of this ensemble member, j of parameter
        end

        % pass namelist to updateParameters function
        Params = updateParameters(Params, FixedParams_s, ensArgIn{:});

    case 'method2' % method 2: load params form hard disk, e.g. by the tag (pSetName)
        if exist(pSetName) 
            loadedParams = load(pSetName);
            Params = loadedParams.Params;
            clear loadedParams
        end
    case 'method3' % method 3: start with a pre-loaded parameter set (like in method 2) and update certain target params
        %load base param set
        pSetName = 'fittedParameters_RMS_Hellinger_ZPratio_Atlantic_correctedDiffusion.mat'

        if exist(pSetName) 
            loadedParams = load(pSetName);
            Params = loadedParams.Params;
            clear loadedParams
        end

        % and update
        switch modTag
            case {'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_fullRunForPaper', ...
                    'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_adjustedVin_fullRunForPaper'}
                newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
                Params = updateParameters(Params, FixedParams_s, 'wPOM1', 10, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2', 'rDON', 0.04, 'rDOC', 0.04)

            case {'optimisedParams_AH_wPOM_100_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_AtlanticOnly_forDiscussion',...
                    'optimisedParams_AH_wPOM_100_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_forDiscussion'}
                % run for discussion, (only atlantic summer), same as
                % reference but higher wPOM
                newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
                Params = updateParameters(Params, FixedParams_s, 'wPOM1', 100, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2', 'rDON', 0.04, 'rDOC', 0.04)

            case 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_1.915e-07_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'
                newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
                Params = updateParameters(Params, FixedParams_s, 'wPOM1', 10, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2', 'rDON', 0.04, 'rDOC', 0.04, 'aP', 1.915e-07)

            case 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_7.66e-07_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'; 
                newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
                Params = updateParameters(Params, FixedParams_s, 'wPOM1', 10, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2', 'rDON', 0.04, 'rDOC', 0.04, 'aP', 7.66e-07)

            case 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_30e-07_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'
                newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
                Params = updateParameters(Params, FixedParams_s, 'wPOM1', 10, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2', 'rDON', 0.04, 'rDOC', 0.04, 'aP', 30e-07)

            case 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_aP_3.83e-06_icelessArcSummer_adjustedVin_ArcticOnly_forDiscussion'
                newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
                Params = updateParameters(Params, FixedParams_s, 'wPOM1', 10, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2', 'rDON', 0.04, 'rDOC', 0.04, 'aP', 3.83e-06)

%             case 'optimisedParams_AH_m_b_-0,1'
%                 % for 'optimisedParams_AH_m_b_-0,1' (decrease m_b)
%                 Params = updateParameters(Params, FixedParams_s, 'm_b', -0.1)
%             case 'optimisedParams_AH_m2_5e-4'
%                 % for 'optimisedParams_AH_m2_5e-4'
%                 Params = updateParameters(Params, FixedParams_s, 'm2', 5e-4)
%             case 'optimisedParams_AH_m2_5e-2'
%                 % for 'optimisedParams_AH_m2_5e-2'
%                 Params = updateParameters(Params, FixedParams_s, 'm2', 5e-2)
%             case 'optimisedParams_AH_m2_5e-2_atLargeSC';
%                 newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2]
%                 Params = updateParameters(Params, FixedParams_s, 'm2', newm2')
%             case 'optimisedParams_AH_m2_0.1_atLargeSC'
%                 newm2 = [0 0 0 0 0 0 0 0 0.1 0 0 0 0 0 0 0 0 0.1]
%                 Params = updateParameters(Params, FixedParams_s, 'm2', newm2')
%             case 'optimisedParams_AH_Qmin_QC_a_0,05'
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_a', 0.05)
%             case 'optimisedParams_AH_Qmin_QC_a_0,01'
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_a', 0.01)
%             case 'optimisedParams_AH_pmax_b_0,01'
%                 % for 'optimisedParams_AH_pmax_b_0,01';
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', 0.01)
%             case 'optimisedParams_AH_pmax_b_-0,01'
%                 % use Aidans opt params, but larger, positive pmax_b:
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.01)
%             case 'optimisedParams_AH_pmax_b_-0,08'
%                 % use Aidans opt params, but larger, positive pmax_b:
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08)
%                 case 'optimisedParams_AH_pmax_b_-0,1'
%                 % use Aidans opt params, but smaller, negative pmax_b:
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.1)
%             case 'optimisedParams_AH_pmax_b_-0,2'
%                 % use Aidans opt params, but smaller, negative pmax_b:
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.2)
%             case 'optimisedParams_AH_Qmin_QC_b_-0,05'
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_b', -0.05)
%             case 'optimisedParams_AH_Qmin_QC_a_0.05';
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_a', 0.05)
%             case 'optimisedParams_AH_Qmin_QC_a_0.05_theta_2';
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_a', 0.05, 'theta', 2)
%             case 'optimisedParams_AH_Qmin_QC_a_0.05_theta_6';
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_a', 0.05, 'theta', 6)
%                 case 'optimisedParams_AH_Qmin_QC_a_0.05_theta_12';
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_a', 0.05, 'theta', 12)
%             case 'optimisedParams_AH_rPOM_0.02'
%                 Params = updateParameters(Params, FixedParams_s, 'rPOC', 0.02, 'rPON', 0.02)
%             case 'optimisedParams_AH_pmax_b_-0,01_Qmin_QC_a_0.05'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.01, 'Qmin_QC_a', 0.05)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05';
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05)
%             case 'optimisedParams_AH_pmax_b_-0,01_Qmin_QC_a_0.05_theta_1,6'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.01, 'Qmin_QC_a', 0.05, 'theta', 1.6)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_1,6'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 1.6)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_3'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 3)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_3_m2_5e-2_atLargeSC'
%                 newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2]
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 3, 'm2', newm2')
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_6_m2_5e-2_atLargeSC'
%                 newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2]
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 6, 'm2', newm2')
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_9_m2_5e-2_atLargeSC'
%                 newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2]
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 9, 'm2', newm2')
%             case 'optimisedParams_AH_pmax_b_-0,01_Qmin_QC_a_0.05_theta_3_m2_5e-2_atLargeSC'
%                 newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2]
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.01, 'Qmin_QC_a', 0.05, 'theta', 3, 'm2', newm2')
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_Qmax_delQ_a_1.4284_Qmax_delQ_b_0'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05) % first update pmax and Qmin_QC
%                 Params = updateParameters(Params, FixedParams_s, 'Qmax_delQ_a', 1.4284, 'Qmax_delQ_b', 0) % then update Qmax_QC
%                 %Qmax is 0.0350, when Qmin is 0.05
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_Qmax_delQ_a_0.2941_Qmax_delQ_b_0'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05) % first update pmax and Qmin_QC
%                 Params = updateParameters(Params, FixedParams_s, 'Qmax_delQ_a', 0.2941, 'Qmax_delQ_b', 0) % then update Qmax_QC
%                 % Qmax is 0.17, when Qmin i 0.05
%             case 'optimisedParams_AH_Qmin_QC_a_0.05_Qmax_delQ_a_0.2941_Qmax_delQ_b_0'
%                 Params = updateParameters(Params, FixedParams_s, 'Qmin_QC_a', 0.05) % first update Qmin_QC
%                 Params = updateParameters(Params, FixedParams_s, 'Qmax_delQ_a', 0.2941, 'Qmax_delQ_b', 0) % then update Qmax_QC
%                 % Qmax is 0.17, when Qmin i 0.05
%             case 'optimisedParams_AH_aP_2e-7'
%                 Params = updateParameters(Params, FixedParams_s, 'aP', 2e-7)
%             case 'optimisedParams_AH_aP_8e-7'
%                 Params = updateParameters(Params, FixedParams_s, 'aP', 8e-7)
%             case 'optimisedParams_AH_aP_8e-6'
%                 Params = updateParameters(Params, FixedParams_s, 'aP', 8e-6)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_6_m2_5e-2_atLargeSC_aP_8e-7'
%                 newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2]
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 6, 'm2', newm2', 'aP', 8e-7)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_aP_6e-7'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'aP', 6e-7)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_aP_6e-7_ArcSummer-temp-ArcAutumn'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'aP', 6e-7)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_aP_6e-7_icelessArcSummer'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'aP', 6e-7)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_aP_6e-7_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'aP', 6e-7)
%              case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
%                  % this is the latest reference solution! 
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6)
            % case 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_icelessArcSummer_adjustedVin_fullRunForPaper'
            %     newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
            %     Params = updateParameters(Params, FixedParams_s, 'wPOM1', 10, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2')
            
            % case 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_rPOM_0.02_icelessArcSummer_adjustedVin'
            %     newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2];
            %     Params = updateParameters(Params, FixedParams_s, 'wPOM1', 10, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11, 'm2', newm2', 'rDON', 0.04, 'rDOC', 0.04, 'rPON', 0.02, 'rPOC', 0.02)

                
%               case 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6)
%             case 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_Gmax_a_11_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'Gmax_a', 11)
%             case 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.04_theta_5.6_Gmax_a_11_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'pmax_b', -0.08, 'Qmin_QC_a', 0.04, 'theta', 5.6, 'Gmax_a', 11)
%             case 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.035_theta_5.6_Gmax_a_11_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 5.6, 'Gmax_a', 11)
            % case 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin'
            %     Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11)
            % case 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin_fullRunForPaper'
            %     Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_0.5_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 0.5)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin_fullRunForPaper'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_icelessArcSummer_adjustedVin_closed'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 11)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_5_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 5)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_22_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 22)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_16_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.035, 'theta', 4.2, 'Gmax_a', 16)
%             case 'optimisedParams_AH_wPOM_50_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 50, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6)
%             case 'optimisedParams_AH_wPOM_5_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'Qmin_QC_a', 0.05, 'theta', 5.6)
%             case 'optimisedParams_AH_wPOM_5_pmax_a_3.5_pmax_b_-0.1_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'pmax_a', 3.5, 'pmax_b', -0.1)
%             case 'optimisedParams_AH_wPOM_5_pmax_b_-0,08_Qmin_QC_a_0.04_theta_5.6_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'wPOM1', 5, 'pmax_b', -0.08, 'Qmin_QC_a', 0.04, 'theta', 5.6)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_icelessArcSummer_adjustedVin_fullRunForPaper'
%                  % this was the latest reference solution before i fixed WPOM1! 
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6)               
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_aP_5e-7_icelessArcSummer_adjustedVin'
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'aP', 5e-7)
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_m2_5e-2_atLargeSC_icelessArcSummer_adjustedVin'
%                 newm2 = [0 0 0 0 0 0 0 0 5e-2 0 0 0 0 0 0 0 0 5e-2]
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'm2', newm2')
%             case 'optimisedParams_AH_pmax_b_-0,08_Qmin_QC_a_0.05_theta_5.6_m2_5e-1_atLargeSC_icelessArcSummer_adjustedVin'
%                 newm2 = [0 0 0 0 0 0 0 0 5e-1 0 0 0 0 0 0 0 0 5e-1]
%                 Params = updateParameters(Params, FixedParams_s, 'pmax_b', -0.08, 'Qmin_QC_a', 0.05, 'theta', 5.6, 'm2', newm2')
        end
end

Params

%% manipulations of forcing
% change aspects of forcing or initial values... for testing


% if modTag contains 'ArcSummer-temp-ArcAutumn', use temperature profiles
% from autumn instead (in Arctic) 
if contains(modTag, 'ArcSummer-temp-ArcAutumn') 
    
    wmfilter_s = strcmp(Forc_s.waterMass, 'Arctic')
   
    
    figure
    p = pcolor(1:304, FixedParams_s.z, Forc_s.T(:,:,wmfilter_s))
    set(p, 'EdgeColor', 'none'); 
        c = colorbar(); 
        c.Label.String = 'Temp degC';
        colormap(jet(30))
         c.YTick = [-2:9]
    %    caxis([-2 9])
    
  
  wmfilter_a = strcmp(Forc_a.waterMass, 'Arctic')
  len_s = size(Forc_s.T,2)
  
    % reassign autumn temps to arctic summer trajs for testing the
    % influence of temperature
  Forc_s.T(:,:,wmfilter_s) = Forc_a.T(:,1:len_s,wmfilter_a)
  
  
   figure
    p = pcolor(1:304, FixedParams_s.z, Forc_s.T(:,:,wmfilter_s))
    set(p, 'EdgeColor', 'none'); 
        c = colorbar(); 
        c.Label.String = 'Temp degC';
        colormap(jet(30))
         c.YTick = [-2:9]
    
  clear wmfilter* len_s
  
end


if contains(modTag, '_closed')  % close system for mass balance testing
    FixedParams_a.POM_is_lost = false;
    FixedParams_s.POM_is_lost = false;
    
    Params.wPOM(end) = 0;
    Params.wk(end,2) = 0;
    Params.wk(end,4) = 0;
    
    Forc_s.K(end,:,:) = 0;
    Forc_a.K(end,:,:) = 0;
    
    
end


%% run model
    % Input data (forcing trajectories and fitting data) may be filtered by the
    % origin of particle trajectories (Atlantic or Arctic)
    % [Forc, Data] = filterInputByOrigin(Forc, Data, 'fitTrajectories', 'Atlantic');

    % Integrate
    % Initialise state variables.
    % Store state variables in array v0: 1st dimension = variable
    %                                    2nd dimension = location (trajectory)
    v0_s = initialiseVariables(FixedParams_s, Params, Forc_s);
    v0_a = initialiseVariables(FixedParams_a, Params, Forc_a);
    
if contains(modTag, 'adjustedVin') 
    % double the DIN values for arctic trajectories in summer and autumn
    % to correct too low SINMOD estimations for the Arctic water
    
    % summer
    v0_s(FixedParams_s.IN_index, strcmp(Forc_s.waterMass, 'Arctic')) = ...
        repmat([10; 10; 10; 10; 10; 10; 11; 11; 11], 1, sum(strcmp(Forc_s.waterMass, 'Arctic')));        %   v0_s(FixedParams_s.IN_index, strcmp(Forc_s.waterMass, 'Arctic')).*2
    % autumn
    v0_a(FixedParams_a.IN_index, strcmp(Forc_a.waterMass, 'Arctic')) = ...
        repmat([11; 11; 11; 11; 11; 11; 11; 11; 11],1, sum(strcmp(Forc_a.waterMass, 'Arctic')));
end   
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Order of variables = inorganic nutrients [depth]
    %                      phytoplankton       [size, depth, nutrient]
    %                      zooplankton         [size, depth, nutrient]
    %                      organic matter      [type, depth, nutrient]
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Parallelise integrations over trajectories
    poolObj = gcp('nocreate');
    if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

    % Run the model
    tic; disp('.. started integration at'); disp(datetime('now'))
    [out_s, auxVars_s] = integrateTrajectories(FixedParams_s, Params, Forc_s, v0_s, ... 
        FixedParams_s.odeIntegrator, FixedParams_s.odeOptions);
 
    [out_a, auxVars_a] = integrateTrajectories(FixedParams_a, Params, Forc_a, v0_a, ... 
        FixedParams_a.odeIntegrator, FixedParams_a.odeOptions);
    toc
    
    
%     % get modData, matching size data by water mass origin
%     modData_s = matchModOutput2Data2(out_s, auxVars_s, Data_s, FixedParams_s, ...
%         'fitToFullSizeSpectra', 'trueByWaterMass', 'Forc', Forc_s, 'evTrajIX', summerTrajsIndex);
%     modData_a = matchModOutput2Data2(out_a, auxVars_a, Data_a, FixedParams_a, ...
%         'fitToFullSizeSpectra', 'trueByWaterMass', 'Forc', Forc_a, 'evTrajIX', autumnTrajsIndex);
%    
    % get modData, matching size data by station, analogously to scalar
    % data
    modData_s = matchModOutput2Data2(out_s, auxVars_s, Data_s, FixedParams_s, ...
        'fitToFullSizeSpectra', true, 'Forc', Forc_s, 'evTrajIX' , summerTrajsIndex);  %
    modData_a = matchModOutput2Data2(out_a, auxVars_a, Data_a, FixedParams_a, ...
        'fitToFullSizeSpectra', true, 'Forc', Forc_a);  % , 'evTrajIX', autumnTrajsIndex

%% save output    

%summer
filename = [Directories.resultsDir  'summer_output_' modTag '.mat'];
    m = matfile(filename, 'Writable', true);
    m.modData = modData_s;
    m.out = out_s; 
    m.auxVars = auxVars_s;
    m.v0 = v0_s;
    m.Forc = Forc_s;
    m.Data = Data_s;
    m.FixedParams = FixedParams_s;
    m.Params = Params; 
    
    if(exist('summerTrajsIndex', 'var'))
        m.summerTrajsIndex = summerTrajsIndex;
    end

    
    % autumn
filename = [Directories.resultsDir  'autumn_output_' modTag '.mat'];
    m = matfile(filename, 'Writable', true);
    m.modData = modData_a;
    m.out = out_a; 
    m.auxVars = auxVars_a;
    m.v0 = v0_a;
    m.Forc = Forc_a;
    m.Data = Data_a;
    m.FixedParams = FixedParams_a;
    m.Params = Params; 


    
% %% some (deprecated) test plots
% 
% % phytoplankton biomass
% size(out_a.P)
% temp = out_a.P(:, :, find(FixedParams_a.PP_C_index), :, 2 )
% size(temp)
% temp(1,9,1,360)
% 
% % plot traj 53 c biomass
% data = out_a.P(:,:,FixedParams_a.PP_C_index,:,2); % first traj, only C  [sizeclass/OM type, depth, var (N, C, Chl), day, traj]
% % summerize over size classes
% data = squeeze(sum(data, 1));
% 
% data10 = data; % append copy of last row so pcolor can use depth edges on y axis
% data10(10,:) = data10(9,:);
% 
% figure
% tst = pcolor(yearday(Forc_a.t(:,1)), FixedParams_a.zw, data10) % mmol C /m3
% set(tst, 'EdgeColor', 'none');
% %ylim([-100 0])
% c = colorbar('Ticks', [0.01 0.1 1 10]);
% c.Label.String = 'Phy C (mmol m^{-3})';
% colormap(jet(30))
% %caxis([0.1 10])
% set(gca,'ColorScale','log');
% xlabel('Yearday')
% ylabel('Depth (m)')
% xticks(100:100:FixedParams_a.nt)
% xticklabels(yearday(Forc_a.t(100:100:FixedParams_a.nt,1)))
% 
% 
% % hier gibt es NaNs! also beim PP carbon in traj 53... (bei traj 55 nur in
% % tiefe 9)
% % das ist aber nicht bei allen ensemble members gleich schlimm
% 
% % zooplankton biomass
% temp = out_a.Z(:,:,find(FixedParams_a.ZP_C_index), :, 2 )
% 
% % POC
% size(out_a.OM)
% temp = out_a.OM(FixedParams_a.DOM_index,:,FixedParams_a.OM_C_index,:, 2 )
% 
% % DOC
% temp = out_a.OM(FixedParams_a.POM_index,:,FixedParams_a.OM_C_index,:, 2 )
% 
% % bei ID 22 passiert irgendwas an zwischen tag 268 und 269!!!!
%     
%     
%     
% 
%     
%     %% make some plots to show model output
%     % lieber in ein extra skript packen...
%     
% % --- FORCING ---     
% % diffusivity mean of both trajs   
% tiledlayout(1,2)    
% logScale = true;
% ColourBar = true;
% unitDiffusivity = 'diffusivity (m^2 day^{-1})';
% nexttile
% plot_forcing_contour_DepthTime('K', Forc_s, auxVars_s, FixedParams_s, ... 
%     'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
%     'logScale', logScale);  
% title('Diffusivity in Summer set-up');
% 
% nexttile
% plot_forcing_contour_DepthTime('K', Forc_a, auxVars_a, FixedParams_a, ... 
%     'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
%     'logScale', logScale);  
% title('Diffusivity in Autumn set-Up');
% linkaxes();
% 
% % diffusivity for both trajs individually   
% tiledlayout(2,2)    
% logScale = true;
% ColourBar = true;
% unitDiffusivity = 'diffusivity (m^2 day^{-1})';
% nexttile
% plot_forcing_contour_DepthTime('K', Forc_s, auxVars_s, FixedParams_s, ... 
%     'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
%     'logScale', logScale, 'traj', 1);  
% title('Summer set-up, atlantic trajectory');
% nexttile
% plot_forcing_contour_DepthTime('K', Forc_s, auxVars_s, FixedParams_s, ... 
%     'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
%     'logScale', logScale, 'traj', 2);  
% title('Summer set-up, arctic trajectory');
% 
% 
% nexttile
% plot_forcing_contour_DepthTime('K', Forc_a, auxVars_a, FixedParams_a, ... 
%     'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
%     'logScale', logScale, 'traj', 1);  
% title('Autumn set-Up, atlantic trajectory');
% nexttile
% plot_forcing_contour_DepthTime('K', Forc_a, auxVars_a, FixedParams_a, ... 
%     'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
%     'logScale', logScale, 'traj', 2);  
% title('Autumn set-Up, arctic trajectory');
% linkaxes();
% 
% 
% 
% 
% 
% %%
% 
% % --- OUTPUT ---
% colourMap = 'plasma';
% % colourMap = 'viridis';
% % colourMap = 'magma';
% % colourMap = 'inferno';
% % colourMap = 'parula';
% 
% % DIN
% ColourBarLabel = 'DIN (mmol N m^{-3})';
% 
% tiledlayout(2,2)
% nexttile
% Title = ['Summer set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('DIN' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% colormap(colourMap);
% nexttile
% Title = ['Summer set-up, arctic trajectory'];
% plot_output_contour_DepthTime('DIN' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% colormap(colourMap);
% nexttile
% Title = ['Autumn set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('DIN' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% colormap(colourMap);
% nexttile
% Title = ['Autumn set-up, arctic trajectory'];
% plot_output_contour_DepthTime('DIN' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% colormap(colourMap);
% linkaxes();
% 
% 
% 
% % POC
% ColourBarLabel = 'POC (mmol C m^{-3})';
% 
% tiledlayout(2,2)
% nexttile
% Title = ['Summer set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('POC' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Summer set-up, arctic trajectory'];
% plot_output_contour_DepthTime('POC' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('POC' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, arctic trajectory'];
% plot_output_contour_DepthTime('POC' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% linkaxes();
% 
% 
% % DOC
% ColourBarLabel = 'DOC (mmol C m^{-3})';
% 
% tiledlayout(2,2)
% nexttile
% Title = ['Summer set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('DOC' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Summer set-up, arctic trajectory'];
% plot_output_contour_DepthTime('DOC' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('DOC' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, arctic trajectory'];
% plot_output_contour_DepthTime('DOC' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% linkaxes();
% 
% % PON
% ColourBarLabel = 'PON (mmol N m^{-3})';
% 
% tiledlayout(2,2)
% nexttile
% Title = ['Summer set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('PON' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Summer set-up, arctic trajectory'];
% plot_output_contour_DepthTime('PON' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('PON' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, arctic trajectory'];
% plot_output_contour_DepthTime('PON' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% linkaxes();
% 
% % POC
% ColourBarLabel = 'DON (mmol N m^{-3})';
% 
% tiledlayout(2,2)
% nexttile
% Title = ['Summer set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('DON' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Summer set-up, arctic trajectory'];
% plot_output_contour_DepthTime('DON' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, atlantic trajectory'];
% plot_output_contour_DepthTime('DON' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Atlantic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% nexttile
% Title = ['Autumn set-up, arctic trajectory'];
% plot_output_contour_DepthTime('DON' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel);
% linkaxes();
% 
% % noch mehr plots
% 
% 
% % Phy N 
% % summer set up, both trajs, per size class
% ColourBarLabel = 'Phy N (mmol N m^{-3})';
% 
% tiledlayout(9,2)
% 
% for sc = 1:FixedParams_s.nPP_size
% nexttile
% plot_output_contour_DepthTime('P_N' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Atlantic', 'ColourBarLabel', ColourBarLabel, ...
%     'sizeClass', sc);
% 
% nexttile
% plot_output_contour_DepthTime('P_N' , out_s , auxVars_s, FixedParams_s, Forc_s, ...
%     'waterOrigin', 'Arctic', 'ColourBarLabel', ColourBarLabel, ...
%     'sizeClass', sc);
% 
% end
% % UPS, die achsenbeschriftung (y, depth) ist verkehrt (labels nicht reversed, werte stimmen aber)
% % !!!!!!!!!!!!!!!!!!!!!!!!
% 
% % -> Aidans functions can plot P and Z output only per size class. 
% % make new ones to plot total P and Z biomass first, then look into size
% % structure
% plot_output_contour_DepthTime('P_N' , out_a , auxVars_a, FixedParams_a, Forc_a, ...
%     'waterOrigin', 'Arctic' , 'Title', Title, 'ColourBarLabel', ColourBarLabel, ...
%     'sizeClass', 9);
% 
% 
% % plot observed and modeled size spectra
% % histogram
% plotModelSizeSpectra(out_s, auxVars_s, FixedParams_s, Data_s, Forc_s, 'histogram',...
%     'abundance', 'phyto', 200, 'averaged', 'waterOrigin', 'Atlantic')
% 
% % spectrum
% % summer
% plotModelSizeSpectra(out_s, auxVars_s, FixedParams_s, Data_s, Forc_s, 'spectrum',...
%     'abundance', 'phyto', 200, 'averaged', 'waterOrigin', 'Atlantic')
% hold on
% ix = strcmp(Data_s.size.scenario, 'S1') & strcmp(Data_s.size.trophicLevel, 'autotroph');
% plot(Data_s.size.ESD(ix), Data_s.size.cellDensity(ix))
% hold off
% 
% %autumn
% plotModelSizeSpectra(out_a, auxVars_a, FixedParams_a, Data_a, Forc_a, 'spectrum',...
%     'abundance', 'phyto', 260, 'averaged', 'waterOrigin', 'Atlantic')
% hold on
% ix = strcmp(Data_a.size.scenario, 'S3') & strcmp(Data_a.size.trophicLevel, 'autotroph');
% plot(Data_a.size.ESD(ix), Data_a.size.cellDensity(ix))
% ylim([1e5, 1e11])
% hold off
% 
% 

