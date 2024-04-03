%% ectract coordinates of trajectories with dates and months

% Forc contains all trajectories

ForcBackup = Forc; 
Forc = Forc0; 

longTable = table; 
for i = 1:Forc.nTraj
    
    % ectract traj coordinates and metadata
    latTraj = Forc.y(:,i); 
    lonTraj = Forc.x(:,i); 
    timeTraj = Forc.t(:,i);
    idTraj = repelem(Forc.iTraj(i), length(latTraj), 1); 
    indexTraj = repelem(i, length(latTraj), 1);
    wmTraj = repelem(Forc.waterMass(i), length(latTraj), 1);
    monthTraj = month(datetime(timeTraj, 'ConvertFrom', 'datenum')); 
    
    % create table
    trajTable = table(latTraj, lonTraj, timeTraj, idTraj, indexTraj, wmTraj, monthTraj); 
    
    % append table to long table of all Trajs
    longTable = [longTable; trajTable]; 
   
end

% save
% save('~/Documents/microARC model/trajCoords', longTable)
% csvwrite('~/Documents/microARC model/trajCoords', longTable)
writetable(longTable, '~/Documents/microARC model/Sat&WOA trajectories/trajCoordsAutumn.csv')

Forc = ForcBackup; 