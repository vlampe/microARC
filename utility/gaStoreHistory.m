function [state,options,optchanged] = gaStoreHistory(options,state,flag)
% Genetic algorithm (ga) output function (OutputFcn) that dynamically
% stores the population history to workspace variable gapopulationhistory.
% This is useful for debugging if errors occur part-way through the ga
% algorithm.
persistent history costHistory
optchanged = false;
switch flag
    case 'init'
        history(:,:,1) = state.Population;
        costHistory(:,1) = state.Score;
        assignin('base','gapopulationhistory',history);
        assignin('base','gacosthistory',costHistory);
    case 'iter'
        % Update the history every x generations.
        x = 1;
        if rem(state.Generation,x) == 0
            ss = size(history,3);
            history(:,:,ss+1) = state.Population;
            costHistory(:,ss+1) = state.Score;
            assignin('base','gapopulationhistory',history);
            assignin('base','gacosthistory',costHistory);
        end
    case 'done'
        % Include the final population in the history.
        ss = size(history,3);
        history(:,:,ss+1) = state.Population;
        costHistory(:,ss+1) = state.Score;
        assignin('base','gapopulationhistory',history);
        assignin('base','gacosthistory',costHistory);
end
