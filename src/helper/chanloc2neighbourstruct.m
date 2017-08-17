function  [channeighbstructmat ept_tfce_nb]= chanloc2neighbourstruct(chanlocs,neighbourdist)


sens = [];
sens.pnt = [[chanlocs.X];[chanlocs.Y];[chanlocs.Z]]';

% function [neighbours,channeighbstructmat]=compneighbstructfromgradelec(sens,neighbourdist)

nsensors = size(sens.pnt,1);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
  dist(i,:) = sqrt(sum((sens.pnt(1:nsensors,:) - repmat(sens.pnt(i,:), nsensors, 1)).^2,2))';
end

% find the neighbouring electrodes based on distance
% later we have to restrict the neighbouring electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);
% electrode istelf is not a neighbour
channeighbstructmat = (channeighbstructmat .* ~eye(nsensors));


% for the etp_tfce toolbox we need a different format
[x,y] = find(channeighbstructmat);
nb = zeros(size(channeighbstructmat,1),1+max(arrayfun(@(x)sum(y==x),unique(y))));
for ch = 1:length(channeighbstructmat)
    a = [ch find(channeighbstructmat(ch,:))];
    nb(ch,:) = [a,zeros(1,size(nb,2)-numel(a))];
end
ept_tfce_nb =nb;