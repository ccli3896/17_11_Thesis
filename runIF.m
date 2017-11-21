function APts = runIF(StampsFile,injE,injI)
% Difference between this and runIF2 is that this is SW; 2 is LW.
%
% StampsFile is the name (including .mat) of the TS file. A cell of
% timestamp vectors.
% APts is the cell saved under name StampsFile.

    tiFr = 1000; % Timeframe in ms
    dt = .5; % in ms
    spSp = .0003; % Spontaneous spiking probability at any timestep
    refrTi = 10; % in ms NOT SPECIFIED IN PAPER, but this is about how long an AP is
    numC = 600; 
    numE = 500;
    %injE = 1.5; % injected current to alter excitability
    %injI = .1;
    capac = 1; % capacitance. Compare to NEURON simulations, maybe.
    
    alphj = 1+.3*(rand(numC,1)); % Makes random leak vector [1,1.3]
    
    AVs = cell(numC,2); % 1 holds v, 2 holds syn current
    APts = cell(numC,1); % holds timestamps for action potentials.
    cMat = dlmread('Data/connMat.dat');
    
    % Making all the zero vectors so space doesn't have to be dynamically
    % allocated
    for ind1 = 1:numC
        for ind2 = 1:2
            AVs{ind1,ind2} = zeros(1,(tiFr/dt)+1);
        end
    end
    
    refr = zeros(1,numC); 
    apFlag = 0;
    
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
                apFlag = 1;
            end
        % For random firing with probability spSp
            if (refr(c)<=0)
                fire = rand(1);
                if (fire<=spSp)
                    refr(c) = refrTi;
                    AVs{c,1}(t) = 0;
                    apTime = t*dt;
                    APts{c} = [APts{c} apTime]; 
                    apFlag = 1;
                end
            end
            
        % Making current if an AP fired; sums past+present
            if (apFlag)
                AVs{c,2} = max([AVs{c,2};makeCurr(apTime,tiFr,0,dt)]);
                apFlag = 0;
            end
            
        % Normal summing over all cells.
            if (refr(c)<=0)
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
                        wjk = -10;%-2
                    else % I to I;
                        wjk = -2; %-10
                    end
                    synSum = synSum + wjk*AVs{c2,2}(t);
                end
            end
            % Adjusting for I vs E cell
            if (c<=numE)
                AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),injE,alphj(c),synSum,dt,capac);
            else
                AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),injI,alphj(c),synSum,dt,capac);
            end
            end
            
        end
    end
    save(StampsFile,'APts');
    %pickCell(APts,StampsFile);
    plotCell(APts);
end

function newV = cVolt(prevV,inj,alphj,synSum,dt,capac)
% Calculates the voltage for the next time step.
    % Inputs:
    % prevV is the t-1 voltage for that cell.
    % inj is the injected current for i/e.
    % alphj is pretty much the conductance. Units are 1/R, at least.
    % synSum is the already-calculated summation term.
    % dt and capac are defined and explained in the main function.
    
    coef = prevV-(inj/alphj)-((1/alphj)*synSum);
    newV = (inj/alphj)+(1/alphj)*synSum+coef*exp((-1*dt)/(capac/alphj));
end

function synVec = makeCurr(APtime,tlength,tstart,dt)
    % APtime is the time a spike fires
    % tlength is the timeframe of the whole simulation
    % tstart is 0: the beginning of the entire current vector
    % dt is dt
    %
    % synVec is a vector with the current trace from Jablonski et al. 2007
    
    tvec = linspace(tstart,tlength,((tlength-tstart)/dt)+1);
    synVec = zeros(1,length(tvec));
    
    synVec = exp((-1*(tvec-APtime))/3) - exp((-1*(tvec-APtime))/.3);
    synVec(tvec<APtime) = 0;
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