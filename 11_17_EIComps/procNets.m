function [mns,sds] = procNets(StampsFile)
	% Takes a .mat file of timestamps and returns a vector (numE,1) of means and standard deviations
	% of interspike intervals.

	NUM_E = 500;

	APts = load(StampsFile);
	% Because it loads as a struct:
	APts = APts.APts;

	mns = zeros(NUM_E,1);
	sds = zeros(NUM_E,1);

	for ind = 1:NUM_E
		if (length(APts{ind})>1)
			ISIVec = [];
			for l_i = 1:length(APts{ind})-1
				ISIVec = [ISIVec; APts{ind}(l_i+1)-APts{ind}(l_i)];
			end
			mns(ind) = mean(ISIVec);
			sds(ind) = std(ISIVec);
		end
	end
end