function simil = findSim(TSCell)
% Takes a cell of timestamps, assumed 600 total with 100 I, 500 E, and
% returns a single number representing the similarity of activation across
% clusters. The larger the number, the more similarly the clusters behaved.
% Small numbers mean heterogeneity (good).
%
% Input
% TSCell is a single cell from runIF.
%
% Output
% simil is the returned "combination dot product" value. See journal 11.15.

    NumE = 500;
    NumClust = 5;
    binWidth = 50; % Units of timesteps

   	% Finding the length of the simulation
   	lastTs = zeros(NumE,1);
   	for t_i = 1:NumE
        if(~isempty(TSCell{t_i}))
            lastTs(t_i) = TSCell{t_i}(end);
        end
   	end
   	simTime = max(lastTs);

   	% Counting APs for each bin: constructing the main matrix
   	APCounts = zeros(NumClust,ceil(simTime/binWidth));
   	[NumClust,APcc] = size(APCounts);

   	for c_i = 1:NumE
   		clustI = ceil(c_i/(NumE/NumClust)); % Finds the cluster number by dividing by cluster size
   		for a_i = 1:length(TSCell{c_i})
   			tbinI = ceil(TSCell{c_i}(a_i)/binWidth);
   			APCounts(clustI,tbinI) = APCounts(clustI,tbinI) + 1;
   		end
   	end

   	% Normalizing each COLUMN and calculating the sort-of dot product simultaneously
    simil = 0;
   	for t_i = 1:APcc
   		APCounts(:,t_i) = APCounts(:,t_i)./max(APCounts(:,t_i));
   		dProd = 1;
        if (~isnan(sum(APCounts(:,t_i))))
            for c_i = 1:NumClust
                dProd = dProd*APCounts(c_i,t_i);
            end
        end
   		simil = simil + dProd;
   	end
end
