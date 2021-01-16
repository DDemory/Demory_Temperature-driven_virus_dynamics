%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% Figure 5 -- Ecological states biogeogrpahy given SST projection at global scale.
% David Demory -- Jan 2021

%% SST projections
% Download Temperature % https://cera-www.dkrz.de/WDCC/ui/cerasearch/q?query=*:*&page=0&hierarchy_steps_ss=GFDL_CM2.1_SRESA2_1&entry_type_s=dataset
% download: http://www.ipcc-data.org/
% http://www.ipcc-data.org/sim/gcm_monthly/SRES_AR4/index.html

file_name = 'data/GFCM21_SRA2_1_N_ts.1-1200.nc';
%get info:
%ncdisp(file_name)

% get time (years)
timep = ncread(file_name,'time');

% get SST in celcius
sstp = ncread(file_name,'ts')-273.15;

% get dates
dtp = datetime(timep*24*3600, 'ConvertFrom', 'epochtime', 'Epoch', '1960-01-01');

%% Temperature present (2010 - 2020)
% T present
% first and last months of the present period
id1p = find(dtp == '03-Jan-2010 12:00:00');
id2p = find(dtp == '01-Jan-2020 12:00:00');
% average temperature for the present period
sstmeanp = nanmean(sstp(:,:,id1p:id2p),3); % temp mean from 2010 to 2020;
% clean the data (numerical errors and land locations)
sstmeanp(sstmeanp <= -1.8) = -1.8;
sstmeanp(sstmeanp == -Inf) = NaN;
% get longitudes and latitudes
lonp = ncread(file_name,'lon');
latp = ncread(file_name,'lat');
% SST data for the present period (in double)
SSTglobp = double(sstmeanp);
% grid data and clean-up (important for the m_map package)
[LG,LT]=meshgrid(0:360,-89:89);
mapP=ones(size(LG));
mapP = griddata(lonp,latp,SSTglobp',LG,LT,'cubic');
ind=[1:331 1:60]; % Move left side to right
mapP=mapP(:,ind);
LT=LT(:,ind);
LG=LG(:,ind);
LG(LG>1)=LG(LG>1)-360; %...and subtract 360 to some longitudes

%% Temperature future (2090 - 2020)
% future
% first and last months of the future period
id1f = find(dtp == '13-Jan-2090 00:00:00');
id2f = find(dtp == '11-Nov-2100 12:00:00');
% average temperature for the future period
sstmeanf = nanmean(sstp(:,:,id1f:id2f),3);
% clean the data (numerical errors and land locations)
sstmeanf(sstmeanf <= -1.8) = -1.8;
sstmeanf(sstmeanf == -Inf) = NaN;
% SST data for future period
SSTglobf = double(sstmeanf);
% grid data and clean-up (important for the m_map package)
[LG,LT]=meshgrid(0:360,-89:89);
mapF=ones(size(LG));
mapF = griddata(lonp,latp,SSTglobf',LG,LT,'cubic');
ind=[1:331 1:60]; % Move left side to right
mapF=mapF(:,ind);
LT=LT(:,ind);
LG=LG(:,ind);
LG(LG>1)=LG(LG>1)-360;

%% Plots
% Thermal niche data (pairs B,A and C)
Topt_s = [25.30,25.328,25.382]; % Optimal temperatures of growth
Tmax_s = [33.46,33.474,33.497]; % Maximal temperature of growth
TR0_s  = [28.59,23.63,25.741];  % Optimal R0 temperature
% color
colExt = [255 102 102]/255;     % Habitat loss
colRfg = [60,179,113]/255;      % Refuge
colEpi = [135,206,250]/255;     % Endemic

hfig = figure('color','white','position',[0 0 2000 2000]);

% Subplot locations
cp = [1,7,13];  % present maps
cf = [3,9,15];  % future maps
cb = [5,11,17]; % bar plots
cd = [6,12,18]; % dynamics
H = [1,1,1];    % Habitat loss initiation
R = [12,6,8];   % Refuge initiaiton
E = [13,7,9];   % Endemic initiation

% Loop over the pairs
for s = 1:3;
    
    % Evaluate the states
    Topti = Topt_s(s);
    TR0   = TR0_s(s);
    Tmax  = Tmax_s(s);
    
    % Estimate the ecological state matrix for present period
    ecoTp = mapP;
    ecoTp(mapP<=TR0) = 2;  % Endemic stable
    ecoTp(mapP>TR0) = 1;   % Refuge
    ecoTp(mapP>=Tmax) = 0; % Habitat loss
    
    % present maps
    ax = subplot(3,6,cp(s):cp(s)+1);
    % package m_map
    m_proj('robinson','lon',[-358 0]) % projection robinson
    % pcolor
    m_pcolor(LG,LT,ecoTp);
    % world coast
    m_coast('patch',[250,240,230]/255,'edgecolor','k');
    % axis limits
    if s == 3
        m_grid('tickdir','out','LineWidth',2,'FontSize',8,'box','on','xlabeldir','end');
    else
        m_grid('tickdir','out','LineWidth',2,'FontSize',8,'box','on','Xticklabels','');
    end
    % create a colormap for the 3 states
    colormap(ax,[colExt;colRfg;colEpi])
    caxis([0 2])
    % Turn the angles of xtick vertical
    xtickangle(45)
    % titles
    if s == 1
        title({'Present';'2010-2020'},'fontsize',12);
        text(-6.1,0,'Mic-B/MicV-B','FontSize',18,'fontweight','bold')
    elseif s == 2
        text(-6.1,0,'Mic-A/MicV-A','FontSize',18,'fontweight','bold')
    elseif s == 3
        text(-6.1,0,'Mic-C/MicV-C','FontSize',18,'fontweight','bold')
    end
    
    % Zoomin
    %if s == 1
    %    ax_z = axes('position',[0.33 0.85 .15 .2])
    %elseif s == 2
    %    ax_z = axes('position',[0.33 0.55 .15 .2])
    %elseif s == 3
    %    ax_z = axes('position',[0.33 0.25 .15 .2])
    %end
    %m_proj('robinson','lon',[-330 -150],'lat',[0,30])
    %m_pcolor(LG,LT,ecoTp);
    %m_coast('patch',[250,240,230]/255,'edgecolor','k');
    %colormap(ax_z,[colExt;colRfg;colEpi])
    %m_grid('tickdir','out','LineWidth',2,'FontSize',8,'box','on','Xticklabels','','Yticklabels','');
    %caxis([0 2])
    %xtickangle(45)
    %axis tight
    
    
    
    
    %% Future
    % Estimate the ecological state matrix for present period
    ecoTf = mapF;
    ecoTf(mapF<=TR0) = 2;  % Endemic stable
    ecoTf(mapF>TR0) = 1;   % Refuge
    ecoTf(mapF>=Tmax) = 0; % Habitat loss
    
    % Future maps
    ax = subplot(3,6,cf(s):cf(s)+1);
    % package m_map
    m_proj('robinson','lon',[-358 0]) % robinso projection
    % pcolor
    m_pcolor(LG,LT,ecoTf);
    % world coast
    m_coast('patch',[250,240,230]/255,'edgecolor','k');
    % axis limits
    if s == 3
        m_grid('tickdir','out','LineWidth',2,'FontSize',8,'box','on','xlabeldir','end','YtickLabels','');
    else
        m_grid('tickdir','out','LineWidth',2,'FontSize',8,'box','on','Xticklabels','','Yticklabels','');
    end
    % colormap
    colormap(ax,[colExt;colRfg;colEpi])
    caxis([0 2])
    % xtick angle
    xtickangle(45)
    % title
    if s == 1
        title({'Future';'2090-2100'},'fontsize',12);
    end
    
    %% Bar plots
    
    subplot(3,6,cb(s))
    
    % State coverage in present
    MATp = ecoTp(:);
    % total ocean coverage for present period
    np = length(MATp(isnan(MATp)==0));
    % Habitat loss coverage
    n0p = length(MATp(MATp == 0));
    % Refuge loss coverage
    n1p = length(MATp(MATp == 1));
    % Endemic coverage
    n2p = length(MATp(MATp == 2));
    % Coverage percentages
    pn2p = n2p*100/np
    pn1p = n1p*100/np
    pn0p = n0p*100/np
    
    % State repartition in future
    MATf = ecoTf(:);
    % total ocean coverage for future period
    nf = length(MATf(isnan(MATf)==0));
    % Habitat los coverage
    n0f = length(MATf(MATf == 0));
    % Refuge loss coverage
    n1f = length(MATf(MATf == 1));
    % Endemic coverage
    n2f = length(MATf(MATf == 2));
    % Coverage percentages
    pn2f = n2f*100/nf
    pn1f = n1f*100/nf
    pn0f = n0f*100/nf
    
    % barplots
    bp = barh([pn2f-pn2p,pn1f-pn1p,pn0f-pn0p],'EdgeColor','none');
    set(gca,'TickDir','out')
    box off
    bp.FaceColor = 'flat';
    bp.CData(1,:) = colEpi;
    bp.CData(2,:) = colRfg;
    bp.CData(3,:) = colExt;
    yticklabels({'Endemic','Refuge','Habitat loss'})
    xlim([-15,15])
    xticks([-15,-7.5,0,7.5,15])
    xticklabels({' '})
    % Values on the plots
    text(0+1,1,num2str(E(s)))
    text(pn1f-pn1p+1,2,num2str(R(s)))
    text(pn0f-pn0p+1,3,num2str(H(s)))
    
    % titles
    if s == 1
        title({'% of change between'; 'present en future periods'})
    end
    % Xticks
    if s == 3
        xticks([-15,-7.5,0,7.5,15])
        xticklabels({'-15','-7.5','0','7.5','15'})
        xlabel('%')
    end
    
    
    
    %% Coverage Dynamics
    % Load coverage dynamics
    load('Endemic.mat')
    load('HL.mat')
    load('Refuge.mat')
    
    % dynamics plots
    subplot(3,6,cd(s))
    time = 1:1:1200; % from month 1 (2001)  to month 1200 (2100)
    hold on
    % Endemic dynamic
    plot(Endemic(:,s),'.','MarkerSize',15,'color',colEpi) % model outputs
    mdl_end = fitlm(time,Endemic(:,s)); % Linear regression fit
    plot(time,mdl_end.Coefficients.Estimate(1)+mdl_end.Coefficients.Estimate(2)*time,'LineWidth',2,'color','k')
    ylim([0 100])
    % Refuge dynamic
    plot(Refuge(:,s),'.','MarkerSize',15,'color',colRfg)
    mdl_rfg = fitlm(time,Refuge(:,s));
    plot(time,mdl_rfg.Coefficients.Estimate(1)+mdl_rfg.Coefficients.Estimate(2)*time,'--','LineWidth',2,'color','k')
    % Habitat loss dynamic
    plot(HL(:,s),'.','MarkerSize',15,'color',colExt)
    mdl_hl = fitlm(time,HL(:,s));
    plot(time,mdl_hl.Coefficients.Estimate(1)+mdl_hl.Coefficients.Estimate(2)*time,':','LineWidth',2,'color','k')
    
    % Yticks
    yticks([0 25 50 75 100])
    yticklabels({'0','25','50','75','100'})
    % Xticks
    if s == 3
        xticks([1,290,590,890,1190])
        xticklabels({'05-Jan-2001','29-Jan-2025','23-Jan-2050','17-Jan-2075','11-Jan-2100'});
        xtickangle(45)
        xlabel('Time')
    else
        xticklabels({' '})
    end
    
    % title
    if s == 1
        title({'% of ocean';'coverage'})
    end
    
    
end

set(findall(gcf,'-property','FontSize'),'FontSize',12)

% Save figure
print(hfig,'-depsc','-r600','outputs/Figure5_map.eps')


