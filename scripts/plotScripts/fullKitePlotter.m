%% animations plots
parseLogsout
make_video = 0;
LS = lengthScale;

% plotting time interval
resampleDataRate = 4*sqrt(LS);

% % % video setting
video = VideoWriter('vid_Test', 'Motion JPEG AVI');
video.FrameRate = 20*1/resampleDataRate;

signals = fieldnames(tsc);
timeAnim = 0:resampleDataRate:tsc.(signals{1}).Time(end);

for ii = 1:length(signals)
    tscResample.(signals{ii}) = resample(tsc.(signals{ii}),timeAnim);
end
sol_Rcm_o = squeeze(tscResample.inertialCmPos.Data);
sol_Vcmo = squeeze(tscResample.inertialCmVel.Data);
sol_euler = squeeze(tscResample.eulerAngles.Data);

nNodes = thr.numNodes.Value;
nTethers = numTethers;
n_steps = length(timeAnim);

s_R = cell(nTethers,1);
s_Rn_o = cell(nTethers,1);
s_R1_o = cell(nTethers,1);
sol_outline = cell(5,1);
int_mat1 = zeros(5,3);
int_mat2 = zeros(5,3,n_steps);

for ii = 1:nTethers
    if nTethers > 1
        s_R{ii} = squeeze(tscResample.allNodePos.Data(:,:,ii,:));
        s_R1_o{ii} = s_R{ii}(:,1,:);
        s_Rn_o{ii} = s_R{ii}(:,end,:);
    else
        s_R{ii} = squeeze(tscResample.allNodePos.Data(:,:,:));
        s_R1_o{ii} = s_R{ii}(:,1,:);
        s_Rn_o{ii} = s_R{ii}(:,end,:);
    end
    
end

surfNames = fieldnames(vhcl.surfaceOutlines);

rotSeq = NaN(3,3,n_steps);
Ri_cm = NaN(5,3,n_steps);

for ii = 1:n_steps
    rotSeq(:,:,ii) = rotation_sequence(sol_euler(:,ii));
end
for jj = 1:5
    for ii = 1:n_steps
        for kk = 1:5
            Ri_cm(kk,:,ii) = (sol_Rcm_o(:,ii) + rotSeq(:,:,ii)*vhcl.surfaceOutlines.(surfNames{jj}).Value(:,kk) )';
        end
    end
    sol_outline{jj} = Ri_cm;
    
end


bx = zeros(nTethers,2);
by = zeros(nTethers,2);
bz = zeros(nTethers,2);

plotMargin = 5*LS;

for ii = 1:nTethers
    [xmin,xmax] = bounds(squeeze(s_R{ii}(1,:,:)),'all');
    [ymin,ymax] = bounds(squeeze(s_R{ii}(2,:,:)),'all');
    [zmin,zmax] = bounds(squeeze(s_R{ii}(3,:,:)),'all');
    
    bx(ii,:) = [xmin - mod(xmin,plotMargin) - plotMargin,...
        xmax - mod(xmax,plotMargin) + plotMargin];
    by(ii,:) = [ymin - mod(ymin,plotMargin) - plotMargin,...
        ymax - mod(ymax,plotMargin) + plotMargin];
    bz(ii,:) = [zmin - mod(zmin,plotMargin) - plotMargin,...
        zmax - mod(zmax,plotMargin) + plotMargin];
end

%% plot
if exist('fn','var')
    fn = fn+1;
else
    fn = 1;
end
figure(fn)
set(gcf,'Position',[200 100 2*560 2*420]);

% colors
red = 1/255*[228,26,28];
black = 1/255*[0,0,0];
line_wd = 1;

mov(1:n_steps)=struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren');

for ii = 1:n_steps
    
    if ii > 1
        h = findall(gca,'type','line','color',red,'-or','color',black);
        delete(h);
    end
    
    for kk = 1:nTethers
        p3d_1 = plot3(s_R{kk}(1,:,ii),s_R{kk}(2,:,ii),s_R{kk}(3,:,ii),...
            '-+','linewidth',line_wd,'color',black);
        hold on
        pRcm_n = plot3([s_R{kk}(1,end,ii) sol_Rcm_o(1,ii)],...
            [s_R{kk}(2,end,ii) sol_Rcm_o(2,ii)],...
            [s_R{kk}(3,end,ii) sol_Rcm_o(3,ii)],...
            '-','linewidth',line_wd,'color',red);
        for jj = 1:5
            p_kite = plot3(sol_outline{jj}(:,1,ii),...
                sol_outline{jj}(:,2,ii),...
                sol_outline{jj}(:,3,ii),...
                '-','linewidth',line_wd,'color',red);
            p_fuse = plot3([sol_outline{1}(1,1,ii); sol_outline{end}(1,1,ii)],...
                [sol_outline{1}(1,2,ii); sol_outline{end}(1,2,ii)],...
                [sol_outline{1}(1,3,ii); sol_outline{end}(1,3,ii)],...
                '-','linewidth',line_wd,'color',red);
            
        end
    end
    
    if ii == 1
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
        xlim([-max(abs(bx(:)))-(10*LS) max(abs(bx(:)))+(10*LS)]);
        ylim([-max(abs(by(:))) max(abs(by(:)))]);
        zlim([0 max(bz(:)) + (5*LS)]);
%         axis equal
        hold on
        grid on
    end
    
    title(['Time = ',sprintf('%0.2f', timeAnim(ii)),' s'])
    
    try
%         waitforbuttonpress
    catch
        break
    end
    F(ii) = getframe(gcf);
%     im = frame2im(F(ii));
%     [imind,cm] = rgb2ind(im,256);
%     imwrite(imind,cm,'gif_file','gif','WriteMode','append');
    

end


if make_video == 1
    open(video)
    for i = 1:length(F)
        writeVideo(video, F(i));
    end
    close(video)
end