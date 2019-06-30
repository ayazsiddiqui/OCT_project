%% animations plots

nNodes = tp.tethers(1).numNodes;
nTethers = length(tp.tethers);

s_R = cell(nTethers,1);
s_Rn_o = cell(nTethers,1);
s_R1_o = cell(nTethers,1);

for ii = 1:nTethers
    s_R{ii} = squeeze(tscResample.allNodePos.Data(:,end-1:end,ii,:).*(1/Lscale));
    s_R1_o{ii} = s_R{ii}(:,1,:);
    s_Rn_o{ii} = s_R{ii}(:,end,:);
end

bx = zeros(nTethers,2);
by = zeros(nTethers,2);
bz = zeros(nTethers,2);

rf = 10;

for ii = 1:nTethers
    [xmin,xmax] = bounds(squeeze(s_R{ii}(1,:,:)),'all');
    [ymin,ymax] = bounds(squeeze(s_R{ii}(2,:,:)),'all');
    [zmin,zmax] = bounds(squeeze(s_R{ii}(3,:,:)),'all');
    
    bx(ii,:) = [-(abs(xmin)-mod(abs(xmin),rf)+rf), abs(xmax)-mod(abs(xmax),rf)+rf];
    
    by(ii,:) = [-(abs(ymin)-mod(abs(ymin),rf)+rf), abs(ymax)-mod(abs(ymax),rf)+rf];
    
    bz(ii,:) = [(abs(zmin)-mod(abs(zmin),rf)+rf), abs(zmax)-mod(abs(zmax),rf)+rf];
    
end

%% plot
n_steps = length(time);
if plot_animation == 0
%     return
end
fn = fn+1;
figure(fn)

% % % video setting
video = VideoWriter('vid_Test', 'Motion JPEG AVI');
video.FrameRate = 1/resampleDataRate;

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
        pRcm_n = plot3([s_R{kk}(1,end,ii) sol_Rcm_o(1,1,ii)],...
            [s_R{kk}(2,end,ii) sol_Rcm_o(2,1,ii)],...
            [s_R{kk}(3,end,ii) sol_Rcm_o(3,1,ii)],...
            '-','linewidth',line_wd,'color',red);
    end
    
    if ii == 1
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
        xlim([min(bx(:)) max(bx(:))]);
        ylim([min(by(:)) max(by(:))]);
        zlim([min(bz(:)) max(bz(:))]);
%         axis equal
        hold on
        grid on
    end
    
    title(['Time = ',sprintf('%0.2f', time(ii)),' s'])
    
    try
%         waitforbuttonpress
    catch
        break
    end
    F(ii) = getframe(gcf);

end

if make_video == 1
    open(video)
    for i = 1:length(F)
        writeVideo(video, F(i));
    end
    close(video)
end