%%%%%%%%%%%%%%%%%% Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%% Surplus and DSE 2nd Step parameters %%%%%%%%%%%%%%%%%%%%%%
budgetRT = 25900;   % (us) x
budgetLUTs = 9403;  % y
budgetFFs = 35083;  % z
budgetDSPs = 86;    % Not considered in 3d plot
budgetBRAMs = 29    % Not considered in 3d plot
LBConstr = -100;    % For plotting issues

p_i = 13;   % Number of total processors in the selected solution
phi_i = 8;  % Number of total physical links in the selected solution

metrics = 2;        % Number of total metrics (i.e., cache, bandwidth)
constraints = 5;    % Number of total constraints (i.e., rt, LUT, FF, DSP, BRAM)

m = p_i * metrics;  % Number of total elements considered for the metrics
n = 9;              % Number of unique oCMS

%%%%%%%%%%%%%%%% oCMS Query Results %%%%%%%%%%%%%%%%%%%%%%

oCMS_name = {   'Ho et al. (2014)', ...
                'L3SU (2023)', ...
                'Moro et al. (2015)', ...
                'Valente et al. (2021)', ...
                'Valente et al. (2016)', ...
                'AXI PMU (2023)', ...
                'Kyung et al. (2010)', ...
                'Brilli et al. (2022)', ...
                'Brilli at al. (2023)'}     % Name of unique oCMS

% oCMS_i = [pd lut ff dsp bram] % Unique oCMS (i.e., Table 8)
oCMS = [0  886  303 0 0;    % p3, cache
        0   78   47 0 0;    % p3, cache, bw   % OLD: 333 188 0 0
        0  160  320 0 0;    % p3, cache
        0  887 1014 0 2;    % spp, bw
        0  260  205 0 0;    % p3, bw
     2848 1127 1681 0 1;    % \phi5, bw
        0 4815 3168 0 0;    % \phi5, bw
      300    0    0 0 0;    % p3, bw
        1  397  917 0 0];   % p3, \phi5, bw

%%%%%%%%%%%%%%%% oCMS Frequencies Matrices %%%%%%%%%%%%%%%%%%%%%%
M1_freq_name = {'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p1 bw','p2 bw','p3 bw','p4 bw','p5 bw','p6 bw','p7 bw','p8 bw','p9 bw','p10 bw','p11 bw','p12 bw','p13 bw'}

% Processor clock or Physical Link Transmission frequencies
working_freq = [1/75; 1/75; 1/75; 1/75; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200;  % cache
                1/75; 1/75; 1/75; 1/75; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200]; % bandwidth

M1_freq_rt = [ 0 0 0 0 0 0 0 0 0; % p1 (8051) - cache
            0 0 0 0 0 0 0 0 0; % p2 (8051) - cache
            1 1 1 0 0 0 0 0 0; % p3 (LEON3) - cache
            0 0 0 0 0 0 0 0 0; % p4 (LEON3) - cache
            0 0 0 0 0 0 0 0 0; % p5 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p6 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p7 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p8 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p9 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p10 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p11 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p12 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p13 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p1 (8051) - bw
            0 0 0 0 0 0 0 0 0; % p2 (8051) - bw
            0 1 0 0 1 0 0 1 0; % p3 (LEON3) - bw
            0 0 0 0 0 0 0 0 0; % p4 (LEON3) - bw
            0 0 0 1 0 0 0 0 0; % p5 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p6 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p7 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p8 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p9 (SPP) - bw
            0 0 0 0 0 0 0 0 0; % p10 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p11 (SPP) - bw
            0 0 0 0 0 0 0 0 0; % p12 (SPP) - bw
            0 0 0 0 0 0 0 0 0];  % p13 (SPP) - bw

%%%%%%%%% NOTE: If the element does not scale with frequency, set its frequency to 1

M1_freq_area = [ 0 0 0 0 0 0 0 0 0; % p1 (8051) - cache
            0 0 0 0 0 0 0 0 0; % p2 (8051) - cache
            1 1 1 0 0 0 0 0 0; % p3 (LEON3) - cache
            0 0 0 0 0 0 0 0 0; % p4 (LEON3) - cache
            0 0 0 0 0 0 0 0 0; % p5 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p6 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p7 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p8 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p9 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p10 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p11 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p12 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p13 (SPP) - cache
            0 0 0 0 0 0 0 0 0; % p1 (8051) - bw
            0 0 0 0 0 0 0 0 0; % p2 (8051) - bw
            0 1 0 0 1 0 0 1 0; % p3 (LEON3) - bw
            0 0 0 0 0 0 0 0 0; % p4 (LEON3) - bw
            0 0 0 1 0 0 0 0 0; % p5 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p6 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p7 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p8 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p9 (SPP) - bw
            0 0 0 0 0 0 0 0 0; % p10 (SPP) - bw
            0 0 0 1 0 0 0 0 0; % p11 (SPP) - bw
            0 0 0 0 0 0 0 0 0; % p12 (SPP) - bw
            0 0 0 0 0 0 0 0 0];  % p13 (SPP) - bw

%%%%%%%%% NOTE: If the element does not scale with frequency, set its frequency to 1

%%%%%%%%%%%%%%%% M1 Creation  %%%%%%%%%%%%%%%%%%%%%%

M1 = cell(m,n,constraints);

for i = 1:m
    for j = 1:n
        M1_rt_tmp = M1_freq_rt(i,j) * working_freq(i) * oCMS(j,1);
        M1_area_tmp = M1_freq_area(i,j) * oCMS(j,2:end); 
        M1{i,j} = cat(2,M1_rt_tmp,M1_area_tmp);
    end
end

% Find non-zero indices for each row
nonZeroIndices = cell(1, m);
for i = 1:m
    nonZeroIndices{i} = find(M1_freq_rt(i,:));
end

nonZeroIndicesBackup = nonZeroIndices;
indexSol = find(~cellfun(@isempty,nonZeroIndices));
nonZeroIndices = nonZeroIndices(~cellfun(@isempty, nonZeroIndices));

% Generate combinations of indices
indexCombinations = cell(1, size(nonZeroIndices(1,:),2));
[indexCombinations{:}] = ndgrid(nonZeroIndices{:});
combinations = cell2mat(cellfun(@(x) x(:), indexCombinations, 'UniformOutput', false));

archNum = size(combinations,2)          % number of considered not empty structural elements
totalSolutions = size(combinations,1)   % total number of possible solutions

minVal = 1000000;
feasibleSol = 0;
for j = 1:totalSolutions
    pd(j) = 0;
    luts(j) = 0;
    ffs(j) = 0;
    dsps(j) = 0;
    brams(j) = 0;
    for i = 1:archNum
        indx_row = indexSol(1,i);
        indx_column = combinations(j,i);
        sol_tmp = cell2mat(M1(indx_row,indx_column));
        pd(j) = pd(j) + sol_tmp(1);
        luts(j) = luts(j) + sol_tmp(2);
        ffs(j) = ffs(j) + sol_tmp(3);
        dsps(j) =  dsps(j) + sol_tmp(4);    
        brams(j) = brams(j) + sol_tmp(5);   
    end
    if pd(j) < budgetRT & luts(j) < budgetLUTs & ffs(j) < budgetFFs & dsps(j) < budgetDSPs & brams(j) < budgetBRAMs
        feasibleSol = feasibleSol + 1;
    end
    tot = pd(j) + luts(j) + ffs(j) + dsps(j) + brams(j);  % Linear Cost Function
    if tot <= minVal 
        BestMinVal = tot;
        BestIndexComb = j;
        BestPD = pd(j);
        BestLUT = luts(j);
        BestFF = ffs(j);
        BestDSP = dsps(j);
        BestBRAMS = brams(j);
    end
end

%%%%%%%%%%%%%%%%%%%% Plot Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%

coord=[LBConstr  LBConstr    LBConstr;
       budgetRT  LBConstr    LBConstr;
       budgetRT  budgetLUTs  LBConstr;
       LBConstr  budgetLUTs  LBConstr;
       LBConstr  LBConstr    budgetFFs;
       budgetRT  LBConstr    budgetFFs;
       budgetRT  budgetLUTs  budgetFFs;
       LBConstr  budgetLUTs  budgetFFs;];

idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';

xc = coord(:,1);
yc = coord(:,2);
zc = coord(:,3);

figure

alfa = 1000;
patch(xc(idx), yc(idx), zc(idx), 'r', 'facealpha', 0.3);
axis([0 max([max(pd)+alfa,budgetRT+alfa]) 0 max([max(luts)+alfa,budgetLUTs+alfa]) 0 max([max(ffs)+alfa,budgetFFs+alfa])]) % pd = max(max(pd),max(pdMax))+2000, LUTs = max(max(luts),max(lutsMax))+2000,  FF = max(max(ffs),max(ffsMax))+2000
view(3);

grid on
hold on

%%%%%%%%%%%%%%%%%%%% Plot oCMS Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%

scatter3(pd(1,:)',luts(1,:)',ffs(1,:)','*','MarkerEdgeColor','r', 'MarkerFaceColor',[.49 1 .63]);

xlabel('PD (us)')
ylabel('LUTs')
zlabel('FFs')

% set(gca,'XScale','log')
% set(gca,'YScale','log')
% set(gca,'ZScale','log')

legend('Constraints','oCMS Solutions','Location','NorthWest');

hold on

%%%%%%%%%%%%%%%%%%%% Print oCMS Best Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Feasible Solutions: \n')
feasibleSol(1)

fprintf('Best Solutions: \n')
for i = 1:size(indexSol,2)
    fprintf('Architectural Element: ');
    M1_freq_name{1,indexSol(i)}

    indexTmp = combinations(BestIndexComb,i);
    fprintf('Selected oCMS: ');
    oCMS_name{1,indexTmp}
end

fprintf('Selected oCMS PD: ');
BestPD

fprintf('Selected oCMS LUTs: ');
BestLUT

fprintf('Selected oCMS FFs: ');
BestFF

fprintf('Selected oCMS DSPs: ');
BestDSP

fprintf('Selected oCMS BRAMs: ');
BestBRAMS

toc