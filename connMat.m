function connMat
    % Makes a matrix of cell# by cell# (600) as a reference for runIF.m,
    % telling which cells are dependent on which. 
    % Columns depend on rows. (3,2) means 3 feeds 2
    % Periodic boundary conditions added 10.3 (E and I are separate
    % circles)
    
    E_cells = 500;
    I_cells = 100;
    PerClust = 100; %Cells per cluster
    NumClust = 5;

    E_n = 5; % Radii of connections to nearby neurons
    I_n = 1;
    
    cMat = zeros(E_cells+I_cells);
    
    % Connecting all neighboring excitatory cells
    for ind = 1:E_cells
        % If very low, close to the edge
        if (ind<=E_n)
            cMat(ind,1:ind) = 1;
            cMat(ind,E_cells-E_n+ind:E_cells) = 1;
        else
            cMat(ind,ind-E_n:ind) = 1;
        end
        
        % If very high, close to the edge
        if (ind>E_cells-E_n)
            cMat(ind,ind:E_cells) = 1;
            cMat(ind,1:E_n-(E_cells-ind)) = 1;
        else
            cMat(ind,ind:ind+E_n) = 1;
        end
    end
    
    % Connecting all neighboring inhibitory cells 
    for ind = E_cells+1:E_cells+I_cells
        % If very low, close to the edge
        if (ind<=E_cells+I_n)
            cMat(ind,E_cells+1:ind) = 1;
            cMat(ind,I_cells-I_n+ind:E_cells+I_cells) = 1;
        else
            cMat(ind,ind-I_n:ind) = 1;
        end
        
        % If very high, close to the edge
        if (ind>E_cells+I_cells-I_n)
            cMat(ind,ind:E_cells+I_cells) = 1;
            cMat(ind,E_cells+1:E_cells+(E_cells+I_cells-ind)+1) = 1;
        else
            cMat(ind,ind:ind+I_n) = 1;
        end
    end
    
    % Connecting each inhibitory cell to five neighboring excitatory ones
    for ind = 1:I_cells
        cMat(((ind-1)*5)+1:(ind*5),ind+E_cells) = 1;
    end
    
    % Connecting each excitatory cell to ten random inhibitory ones 
    for ind = 1:E_cells
        IMat = E_cells+1:E_cells+I_cells;
        ITen = IMat(randperm(length(IMat),10));
        cMat(ITen,ind) = 1;
    end
    
    % Addition of excitatory random connections
    % For a ratio of 1:5 local:long-range connections. 
    for e_ind = 1:E_cells
        for ind = 1:2
            tmpN = randi([1 E_cells],1,1);
            if (cMat(e_ind,tmpN)==1)
                ind = ind-1;
                continue;
            end
            cMat(e_ind,tmpN) = 1;
        end
    end
    
    % Addition of inhibitory random connections
    % For a ratio of 1:1 local:long-range
    for i_ind = 1:I_cells
        for ind = 1:2
            tmpN = randi([E_cells+1 E_cells+I_cells],1,1);
            if (cMat(i_ind+E_cells,tmpN)==1)
                ind = ind-1;
                continue;
            end
            cMat(i_ind+E_cells,tmpN) = 1;
        end
    end
    
    % Cluster creation
    for clust = 0:NumClust-1
        for e_ind = (clust*PerClust)+1:(clust+1)*PerClust
            for ind = 1:8 % Fully connected is 1:PerClust-1
                tmpN = randi([(clust*PerClust)+1 (clust+1)*PerClust],1,1);
                
                if (cMat(e_ind,tmpN)==1)
                    ind = ind-1;
                    continue;
                end
                
                cMat(e_ind,tmpN) = 1;
            end
        end
    end
    
    % Erasing all diagonals
    for c=1:E_cells+I_cells
        cMat(c,c) = 0;
    end
    
    dlmwrite('connMat.dat',cMat,'delimiter','\t');
end
    

