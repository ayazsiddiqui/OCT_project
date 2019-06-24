% post processing
%% colors and linewidth
red = 1/255*[228,26,28];
black = 1/255*[0,0,0];
line_wd = 1;

% parse the logged data signals
parseLogsout

%% resample data
resampleDataRate = 1;
% filename = 'testAnimated.gif';
signals = fieldnames(tsc);

newTimeVec = 0:resampleDataRate:tsc.(signals{1}).Time(end);

for ii = 1:length(signals)
    tscResample.(signals{ii}) = resample(tsc.(signals{ii}),newTimeVec);
end

time = tscResample.inertialCmPos.Time;

sol_Ri_o = tscResample.allNodePos.Data;
sol_Rcm_o = tscResample.inertialCmPos.Data;

nNodes = class_thr(1).numNodes;
nTethers = length(class_thr);

s_R = cell(nTethers,1);
s_Rn_o = cell(nTethers,1);
s_R1_o = cell(nTethers,1);

for ii = 1:nTethers
    s_R{ii} = squeeze(sol_Ri_o(:,:,ii,:));
    s_R1_o{ii} = s_R{ii}(:,1,:);
    s_Rn_o{ii} = s_R{ii}(:,end,:);
end

bx = zeros(nTethers,2);
by = zeros(nTethers,2);
bz = zeros(nTethers,2);

for ii = 1:nTethers
    [xmin,xmax] = bounds(squeeze(s_R{ii}(1,:,:)),'all');
    [ymin,ymax] = bounds(squeeze(s_R{ii}(2,:,:)),'all');
    [zmin,zmax] = bounds(squeeze(s_R{ii}(3,:,:)),'all');
    
    bx(ii,:) = [xmin,xmax];
    by(ii,:) = [ymin,ymax];
    bz(ii,:) = [zmin,zmax];
end

%% plot
n_steps = length(time);

for ii = 1:n_steps
    
    figure(1)
    if ii > 1
        h = findall(gca,'type','line','color',red,'-or','color',black);
        delete(h);
    end
    
    for kk = 1:nTethers
        p3d_1 = plot3(s_R{kk}(1,:,ii),s_R{kk}(2,:,ii),s_R{kk}(3,:,ii),...
            '-+','linewidth',line_wd,'color',black);
        hold on
        pRcm_n = plot3([s_R{kk}(1,end,ii) sol_Rcm_o(1,1,ii)],...
            [s_R{kk}(2,end,ii) sol_Rcm_o(2,1,ii)],...
            [s_R{kk}(3,end,ii) sol_Rcm_o(3,1,ii)],...
            '-','linewidth',line_wd,'color',red);
    end
    
    if ii == 1
        xlabel('Y (m)'); ylabel('Y (m)'); zlabel('Z (m)')
        xlim([min(bx(:)) max(bx(:))]);
        ylim([min(by(:)) max(by(:))]);
        zlim([min(bz(:)) max(bz(:))]);
        hold on
        grid on
    end
    
    title(['Time = ',sprintf('%0.2f', time(ii)),' s'])
    
    try
        waitforbuttonpress
    catch
        break
    end
    
end

