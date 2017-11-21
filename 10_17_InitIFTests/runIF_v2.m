function APts = runIF_v2
% Version 2 uses the analytic solution, with current correction (finds max
% current contribution)
    tiFr = 1000; % Timeframe in ms
    dt = .5; % in ms
    spSp = .0003; % Spontaneous spiking probability at any timestep
    refrTi = 10; % in ms NOT SPECIFIED IN PAPER, but this is about how long an AP is
    numC = 600; 
    numE = 500;
    injE = .7; % injected current to alter excitability
    injI = 1.1;
    capac = 1; % capacitance. Compare to NEURON simulations, maybe.
    
    alphj = 1+.3*(rand(numC,1)); % Makes random leak vector [1,1.3]
    
    AVs = cell(numC,1); % 1 holds v, 2 holds syn current
    APts = cell(numC,1); % holds timestamps for action potentials.
    cMat = dlmread('connMatLW.dat');
    
    % Making all the zero vectors so space doesn't have to be dynamically
    % allocated
    for ind1 = 1:numC
        AVs{ind1} = zeros(1,(tiFr/dt)+1);
    end
    
    refr = zeros(1,numC); 
    
    for t = 2:tiFr/dt
        fprintf('%f\n',t*dt);
        
    % Next section is to update the voltage vectors.
        for c = 1:numC
            
            % If AP is in progress.
            if (refr(c) > 0)
                refr(c) = refr(c)-dt;
                AVs{c,1}(t) = 0;
            end
            % If AP should be fired.
            if (AVs{c,1}(t-1) >= 1)
                refr(c) = refrTi;
                AVs{c,1}(t) = 0;
                apTime = t*dt;
                APts{c} = [APts{c} apTime];
            end
            % For random firing with probability spSp
            if (refr(c)<=0)
                fire = rand(1);
                if (fire<=spSp)
                    refr(c) = refrTi;
                    AVs{c,1}(t) = 0;
                    apTime = t*dt;
                    APts{c} = [APts{c} apTime]; 
                end
            end
            
            % Normal summing over all cells.
            if (refr(c)<=0)            
                % Adjusting for I vs E cell
                if (c<=numE)
                    AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),t,injE,alphj(c),cMat,dt,capac,numC,numE,APts,c);
                else
                    AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),t,injI,alphj(c),cMat,dt,capac,numC,numE,APts,c);
                end
            end
            
        end
    end

    pickCell(APts,'IF_Stamps.dat');
    plotCell(APts);
end

function newV = cVolt(prevV,t,inj,alphj,cMat,dt,capac,numC,numE,APts,c)
% Calculates the voltage for the next time step. A single number.
    % Inputs:
    % prevV is the t-1 voltage for that cell.
    % t is time, to figure out tk vectors.
    % inj is the injected current for i/e.
    % alphj is pretty much the conductance. Units are 1/R, at least.
    % cMat is the connectivity matrix.
    % dt and capac are defined and explained in the main function.
    % numC and numE are for assigning the correct weights.
    % APts is the cell of AP timestamps. 
    % c is the cell number
    
% This version uses the analytical solution to the differential eqn for
% voltage. 

TAUF = .3; % From Zochowski paper
TAUS = 3;

% First makes the tk vector for each cell based on what timestep the
% simulation is on. Saves two AP times, to fix current jumps.
    tk = zeros(numC,2)-1000;
    for c_i = 1:numC
        if (~isempty(APts{c_i})) % if there's an APts vector for this cell
            if (size(APts{c_i})>1)
                tk(c_i,1) = APts{c_i}(end-1) - (dt*(t-1));
            else
                tk(c_i,1) = APts{c_i}(end) - (dt*(t-1));
            end
            tk(c_i,2) = APts{c_i}(end) - (dt*(t-1));
        end
    end

% Calculates the summation term (see journal 10.16)
    synSum = 0;
    for c2 = 1:numC
        % to skip the same cell
        if (c == c2)
            continue;
        end
        % Check for the projecting cells and set wjk
        if (cMat(c2,c)==1)
            if (c2<=numE && c<=numE) % both E
                wjk = 2;
            elseif (c2<=numE && c>numE) % E to I
                wjk = 4;
            elseif (c2>numE && c<=numE) % I to E
                wjk = -2;
            else % I to I;
                wjk = -10; 
            end
            sumTerm = @(tkTest,TAUS,alphj,capac,dt,TAUF) ((TAUS/(alphj*TAUS-capac))*(exp(tkTest/TAUS))*...
                (exp(-dt/TAUS)-exp(-alphj*dt/capac)) - ...
                (TAUF/(alphj*TAUF-capac))*(exp(tkTest/TAUF))*...
                (exp(-dt/TAUF)-exp(-alphj*dt/capac)));
            synSum = synSum + wjk*max([sumTerm(tk(c2,1),TAUS,alphj,capac,dt,TAUF); sumTerm(tk(c2,2),TAUS,alphj,capac,dt,TAUF)]);
        end
    end

    newV = prevV*exp(-alphj*dt/capac) + inj/alphj - (inj/alphj)*exp(-alphj*dt/capac) + synSum;
end

function pickCell(treeName,filename)
	% Saves the network to a .dat file where the first column is cell ID number
	% and the second is timestamps. Cell ID numbers are sorted and so are
	% times.
	%
	% Input:
	% treeName is the timestamp cell. 

		totalThings = 0;

		for ind = 1:length(treeName)
			totalThings = totalThings + length(treeName{ind});
		end
		finalTree = zeros(totalThings,2);
		
		m_ind = 1;
		for ind = 1:length(treeName)
			for l_ind = 1:length(treeName{ind})
				finalTree(m_ind,1) = ind;
				finalTree(m_ind,2) = treeName{ind}(l_ind);
				m_ind = m_ind+1;
			end
		end
		dlmwrite(filename,finalTree,'delimiter','\t');
end

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