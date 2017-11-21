function plotCell(treeName)
    % For plotting all the timestamps
        for ind = 1:length(treeName)
            x{ind} = ones(length(treeName{ind}),1)+.1*(ind-1);
        end
        for ind = 1:length(treeName)
            plot(treeName{ind},x{ind},'.');
            hold on
        end
    % End of plot
end