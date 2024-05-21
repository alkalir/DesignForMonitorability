%%%%%%%%%%%%%%%%%% Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%% Surplus and DSE 2nd Step parameters %%%%%%%%%%%%%%%%%%%%%%

%%% 120 bpm %%%
budgetRT = 201; % (us) x
budgetLUTs = 18462; % y
budgetFFs = 40155; % z
budgetDSPs = 0;    % Not considered in 3d plot
budgetBRAMs = 3;    % Not considered in 3d plot
LBConstr = -100;    % For plotting issues

%%% 80 bpm %%%
budgetRT2 = 4490; % (us) x
budgetLUTs2 = 18462; % y
budgetFFs2 = 40155; % z
budgetDSPs2 = 0;    % Not considered in 3d plot
budgetBRAMs2 = 3;    % Not considered in 3d plot
LB2Constr = -100;    % For plotting issues

%%% 60 bpm %%%
budgetRT3 = 10949; % (us) x
budgetLUTs3 = 18462; % y
budgetFFs3 = 40155; % z
budgetDSPs3 = 0;    % Not considered in 3d plot
budgetBRAMs3 = 3;    % Not considered in 3d plot
LB3Constr = -100;    % For plotting issues

p_i = 13;   % Number of total processors in the selected solution
phi_i = 8;  % Number of total physical links in the selected solution

metrics = 1;        % Number of total metrics (i.e., cache, bandwidth)
constraints = 5;    % Number of total constraints (i.e., rt, LUT, FF, DSP, BRAM)

m = p_i * metrics;  % Number of total elements considered for the metrics
n = 16;              % Number of unique oCMS

%%%%%%%%%%%%%%%% oCMS Query Results %%%%%%%%%%%%%%%%%%%%%%

oCMS_name = { 'Ho et al. (2014)', ...
            'L3SU (2023)', ...
            'Moro et al. (2015)', ...
            'Valente et al. (2021)', ...
            'Valente et al. (2016)', ...
            'AXI PMU (2023)', ...
            'Kyung et al. (2010)', ...
            'Brilli et al. (2022)', ...
            'Brilli at al. (2023)', ...
            'Choi et al. (2005)', ...
            'Khan et all (2014)', ...
            'Meng et al. (2008)', ...
            'Isci et al. (2003)', ...
            'Wang et al. (2009)', ...
            'Zoni et al. (2018)', ...
            'Najem at al. (2017)'}    % Name of unique oCMS

% oCMS_i = [pd lut ff dsp bram]
oCMS = [0 0 0 0 0; % No oCMS from Query
        0 0 0 0 0; % No oCMS from Query 
        0 0 0 0 0; % No oCMS from Query
        0 0 0 0 0; % No oCMS from Query
        0 0 0 0 0; % No oCMS from Query
        0 0 0 0 0; % No oCMS from Query
        0 0 0 0 0; % No oCMS from Query
        0 0 0 0 0; % No oCMS from Query
        0 0 0 0 0; % No oCMS from Query
        33300 78 47 0 0; % p1, power Vecchio: 260 205 0 0
        80000 78 47 0 0; % p1, power Vecchio: 260 205 0 0
        1800000 78 47 0 0; % p1, power Vecchio: 260 205 0 0
        2200000 78 47 0 0; % p1, power Vecchio: 260 205 0 0
        100000 78 47 0 0; % p1, power Vecchio: 260 205 0 0
        0 98 68 0 0; % p1, spp, power
        0 8 8 0 0]; % p1, spp, power Vecchio: 0 8 8 0 0 (182x120 bpm, 272x80 bpm, 363x60 bpm)

%%%%%%%%%%%%%%%% oCMS frequencies %%%%%%%%%%%%%%%%%%%%%%
M5_freq_name = {'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13'}

% Processor clock or Physical Link Transmission frequencies
working_freq = [1/12; 1/12; 1/75; 1/75; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200; 1/200];  % power

M5_freq_rt = [  0	0	0	0	0	0	0	0	0	1	1	1	1	1	14	14; % p1 (8051) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p2 (8051) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p3 (LEON3) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p4 (LEON3) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p5 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p6 (SPP) - power 
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p7 (SPP) - power 
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p8 (SPP) - power 
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p9 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p10 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p11 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p12 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0]; % p13 (SPP) - power 

%%%%%%%%% NOTE: If the element does not scale with frequency, set its frequency to 1

M5_freq_area = [0	0	0	0	0	0	0	0	0	1	1	1	1	1	14	14; % p1 (8051) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p2 (8051) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p3 (LEON3) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p4 (LEON3) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;  % p5 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p6 (SPP) - power 
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p7 (SPP) - power 
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p8 (SPP) - power 
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p9 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p10 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p11 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1;  % p12 (SPP) - power
                0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0]; % p13 (SPP) - power 

%%%%%%%%% NOTE: If the element does not scale with frequency, set its frequency to 1

%%%%%%%%%%%%%%%% M3 Creation  %%%%%%%%%%%%%%%%%%%%%%
 
M5 = cell(m,n,constraints); % oCMS solutions from querying

for i = 1:m
    for j = 1:n
        M5_rt_tmp = M5_freq_rt(i,j) * working_freq(i) * oCMS(j,1);
        M5_area_tmp = M5_freq_area(i,j) * oCMS(j,2:end); 
        M5{i,j} = cat(2,M5_rt_tmp,M5_area_tmp);
    end
end

% Find non-zero indices for each row
nonZeroIndices = cell(1, m);
for i = 1:m
    nonZeroIndices{i} = find(M5_freq_rt(i,:));
end

nonZeroIndicesBackup = nonZeroIndices;
indexSol = find(~cellfun(@isempty,nonZeroIndices));
nonZeroIndices = nonZeroIndices(~cellfun(@isempty, nonZeroIndices));

% Generate combinations of indices
indexCombinations = cell(1, size(nonZeroIndices(1,:),2));
[indexCombinations{:}] = ndgrid(nonZeroIndices{:});
combinations = cell2mat(cellfun(@(x) x(:), indexCombinations, 'UniformOutput', false));

archNum = size(combinations,2) % number of considered not empty structural elements
totalSolutions = size(combinations,1) % total number of possible solutions


BestMinVal = 100000000;
BestMinValpd = 100000000;
BestMinValff = 100000000;
BestMinVallut = 100000000;
BestMinVallutff = 100000000;
BestMinValpdlutff = 100000000;

feasibleSol120bpm = 0;
feasibleSol80bpm = 0;
feasibleSol60bpm = 0;

for j = 1:totalSolutions
    pd(j) = 0;
    luts(j) = 0;
    ffs(j) = 0;
    dsps(j) = 0;
    brams(j) = 0;
    for i = 1:archNum
        indx_row = indexSol(1,i);
        indx_column = combinations(j,i);
        sol_tmp = cell2mat(M5(indx_row,indx_column));
        pd(j) = pd(j) + sol_tmp(1);
        luts(j) = luts(j) + sol_tmp(2);
        ffs(j) = ffs(j) + sol_tmp(3);
        dsps(j) =  dsps(j) + sol_tmp(4);
        brams(j) = brams(j) + sol_tmp(5);
    end
    if pd(j) < budgetRT & luts(j) < budgetLUTs & ffs(j) < budgetFFs & dsps(j) < budgetDSPs & brams(j) < budgetBRAMs
        feasibleSol120bpm = feasibleSol120bpm + 1;
    end
    if pd(j) < budgetRT2 & luts(j) < budgetLUTs2 & ffs(j) < budgetFFs2 & dsps(j) < budgetDSPs2 & brams(j) < budgetBRAMs2
        feasibleSol80bpm = feasibleSol80bpm + 1;
    end
    if pd(j) < budgetRT3 & luts(j) < budgetLUTs3 & ffs(j) < budgetFFs3 & dsps(j) < budgetDSPs2 & brams(j) < budgetBRAMs2
        feasibleSol60bpm = feasibleSol60bpm + 1;
    end
    tot = pd(j) + luts(j) + ffs(j) + dsps(j) + brams(j);  % Linear Cost Function
    if tot <= BestMinVal 
        BestMinVal = tot;
        BestIndexComb = j;
        BestPD = pd(j);
        BestLUT = luts(j);
        BestFF = ffs(j);
        BestDSP = dsps(j);
        BestBRAMS = brams(j);
    end
    totpd = pd(j);
    if totpd  <= BestMinValpd
        BestMinValpd = totpd;
        BestIndexCombpd = j;
        BestPDpd = pd(j);
        BestLUTpd = luts(j);
        BestFFpd = ffs(j);
        BestDSPpd = dsps(j);
        BestBRAMSpd = brams(j);
    end
    totff = ffs(j);
    if totff  <= BestMinValff
        BestMinValff = totff;
        BestIndexCombff = j;
        BestPDff = pd(j);
        BestLUTff = luts(j);
        BestFFff = ffs(j);
        BestDSPff = dsps(j);
        BestBRAMSff = brams(j);
    end
    totlut = luts(j);
    if totlut  <= BestMinVallut
        BestMinVallut = totlut;
        BestIndexComblut = j;
        BestPDlut = pd(j);
        BestLUTlut = luts(j);
        BestFFlut = ffs(j);
        BestDSPlut = dsps(j);
        BestBRAMSlut = brams(j);
    end
    totlutff = luts(j) + ffs(j);
    if totlutff  <= BestMinVallutff
        totpdlutff = pd(j);
        if totpdlutff  <= BestMinValpdlutff & totpdlutff > 0
            BestMinValpdlutff = totpdlutff;
            BestMinVallutff = totlutff;
            BestIndexComblutff = j;
            BestPDlutff = pd(j);
            BestLUTlutff = luts(j);
            BestFFlutff = ffs(j);
            BestDSPlutff = dsps(j);
            BestBRAMSlutff = brams(j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%% Plot Constraints (120 bpm) %%%%%%%%%%%%%%%%%%%%%%%%%%%

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

patch(xc(idx), yc(idx), zc(idx), 'r', 'facealpha', 0.3);
view(3);

grid on

hold on

%%%%%%%%%%%%%%%%%%%% Plot Constraints (80 bpm) %%%%%%%%%%%%%%%%%%%%%%%%%%%

LB2ConstrPD = budgetRT+1;
LB2ConstrLUT = budgetLUTs2; 
LB2ConstrFF = budgetFFs2;

coord2=[LB2ConstrPD  LB2Constr    LB2Constr;  % a
       budgetRT2  LB2Constr    LB2Constr;     % b
       budgetRT2  budgetLUTs2  LB2Constr;     % c
       LB2ConstrPD  budgetLUTs2  LB2Constr;   % d
       LB2ConstrPD  LB2Constr    budgetFFs2;   % e
       budgetRT2  LB2Constr    budgetFFs2;     % f
       budgetRT2  budgetLUTs2  budgetFFs2;     % g
       LB2ConstrPD  budgetLUTs2  budgetFFs2;]; % h

idx2 = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';

xc2 = coord2(:,1);
yc2 = coord2(:,2);
zc2 = coord2(:,3);

patch(xc2(idx2), yc2(idx2), zc2(idx2), 'b', 'facealpha', 0.3);
view(3);

grid on

hold on

%%%%%%%%%%%%%%%%%%%% Plot Constraints (60 bpm) %%%%%%%%%%%%%%%%%%%%%%%%%%%

LB3ConstrPD = budgetRT2+1;
LB3ConstrLUT = budgetLUTs3; 
LB3ConstrFF = budgetFFs3;

coord2=[LB3ConstrPD  LB3Constr    LB3Constr;  % a
       budgetRT3  LB3Constr    LB3Constr;   % b
       budgetRT3  budgetLUTs3  LB3Constr;   % c
       LB3ConstrPD  budgetLUTs3  LB3Constr;   % d
       LB3ConstrPD  LB3Constr    budgetFFs3;   % e
       budgetRT3  LB3Constr    budgetFFs3;   % f
       budgetRT3  budgetLUTs3  budgetFFs3;   % g
       LB3ConstrPD  budgetLUTs3  budgetFFs3;]; % h

idx2 = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';

xc2 = coord2(:,1);
yc2 = coord2(:,2);
zc2 = coord2(:,3);

alfa = 1000;
patch(xc2(idx2), yc2(idx2), zc2(idx2), 'g', 'facealpha', 0.3);
axis([0 15000 0 max([max(luts)+alfa,budgetLUTs3+alfa]) 0 max([max(ffs)+alfa,budgetFFs3+alfa])]) 
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

legend('Constraints 120 bpm','Constraints 80 bpm','Constraints 60 bpm','oCMS Solutions','Location','NorthWest')

hold on

%%%%%%%%%%%%%%%%%%%% Print oCMS Best Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Feasible Solutions 120 bpm: \n')
feasibleSol120bpm(1)

fprintf('Feasible Solutions 80 bpm: \n')
feasibleSol80bpm(1)

fprintf('Feasible Solutions 60 bpm: \n')
feasibleSol60bpm(1)

%%% All constraints %%%
fprintf('Best Solutions: \n')
for i = 1:size(indexSol,2)
    fprintf('Architectural Element: ');
    M5_freq_name{1,indexSol(i)}

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

%%% C(rt) constraint %%%
fprintf('Best Solutions for C(rt): \n')
for i = 1:size(indexSol,2)
    fprintf('Architectural Element: ');
    M5_freq_name{1,indexSol(i)}

    indexTmp = combinations(BestIndexCombpd,i);
    fprintf('Selected oCMS  for C(rt): ');
    oCMS_name{1,indexTmp}
end

fprintf('Selected oCMS PD for C(rt): ');
BestPDpd

fprintf('Selected oCMS LUTs for C(rt): ');
BestLUTpd

fprintf('Selected oCMS FFs for C(rt): ');
BestFFpd

fprintf('Selected oCMS DSPs for C(rt): ');
BestDSPpd

fprintf('Selected oCMS BRAMs for C(rt): ');
BestBRAMSpd

%%% C(lut) constraint %%%
fprintf('Best Solutions for C(lut): \n')
for i = 1:size(indexSol,2)
    fprintf('Architectural Element: ');
    M5_freq_name{1,indexSol(i)}

    indexTmp = combinations(BestIndexComblut,i);
    fprintf('Selected oCMS for C(lut): ');
    oCMS_name{1,indexTmp}
end

fprintf('Selected oCMS PD for C(lut): ');
BestPDlut

fprintf('Selected oCMS LUTs for C(lut): ');
BestLUTlut

fprintf('Selected oCMS FFs for C(lut): ');
BestFFlut

fprintf('Selected oCMS DSPs for C(lut): ');
BestDSPlut

fprintf('Selected oCMS BRAMs for C(lut): ');
BestBRAMSlut

%%% C(ff) constraints %%%
fprintf('Best Solutions for C(ff): \n')
for i = 1:size(indexSol,2)
    fprintf('Architectural Element: ');
    M5_freq_name{1,indexSol(i)}

    indexTmp = combinations(BestIndexCombff,i);
    fprintf('Selected oCMS for C(ff): ');
    oCMS_name{1,indexTmp}
end

fprintf('Selected oCMS PD for C(ff): ');
BestPDff

fprintf('Selected oCMS LUTs for C(ff): ');
BestLUTff

fprintf('Selected oCMS FFs for C(ff): ');
BestFFff

fprintf('Selected oCMS DSPs for C(ff): ');
BestDSPff

fprintf('Selected oCMS BRAMs for C(ff): ');
BestBRAMSff

%%% C(lut) and C(ff) constraints %%%
fprintf('Best Solutions for C(lut) and C(ff): \n')
for i = 1:size(indexSol,2)
    fprintf('Architectural Element: ');
    M5_freq_name{1,indexSol(i)}

    indexTmp = combinations(BestIndexComblutff,i);
    fprintf('Selected oCMS for C(lut) and C(ff): ');
    oCMS_name{1,indexTmp}
end

fprintf('Selected oCMS PD for C(lut) and C(ff): ');
BestPDlutff

fprintf('Selected oCMS LUTs for C(lut) and C(ff): ');
BestLUTlutff

fprintf('Selected oCMS FFs for C(lut) and C(ff): ');
BestFFlutff

fprintf('Selected oCMS DSPs for C(lut) and C(ff): ');
BestDSPlutff

fprintf('Selected oCMS BRAMs for C(lut) and C(ff): ');
BestBRAMSlutff

toc
