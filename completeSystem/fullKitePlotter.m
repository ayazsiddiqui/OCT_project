%% animations plots

nNodes = tp.tethers(1).numNodes;
nTethers = length(tp.tethers);
dsgnData = tp.aeroDesignData;
n_steps = length(time);



s_R = cell(nTethers,1);
s_Rn_o = cell(nTethers,1);
s_R1_o = cell(nTethers,1);
sol_outline = cell(5,1);
int_mat1 = zeros(5,3);
int_mat2 = zeros(5,3,n_steps);

for ii = 1:nTethers
    if nTethers > 1
        s_R{ii} = squeeze(tscResample.allNodePos.Data(:,:,ii,:).*(1/Lscale));
        s_R1_o{ii} = s_R{ii}(:,1,:);
        s_Rn_o{ii} = s_R{ii}(:,end,:);
    else
        s_R{ii} = squeeze(tscResample.allNodePos.Data(:,:,:).*(1/Lscale));
        s_R1_o{ii} = s_R{ii}(:,1,:);
        s_Rn_o{ii} = s_R{ii}(:,end,:);
    end
    
end


for jj = 1:5
    
    for ii = 1:n_steps
        
        for kk = 1:5
            int_mat1(kk,:) = ( sol_Rcm_o(:,ii) + ...
                rotation_sequence(sol_euler(:,ii))...
                *dsgnData.outlines(jj).pts(kk,:)' );
        end
        int_mat2(:,:,ii) = int_mat1;
    end
    sol_outline{jj} = int_mat2;
    
end

bx = zeros(nTethers,2);
by = zeros(nTethers,2);
bz = zeros(nTethers,2);

for ii = 1:nTethers
    [xmin,xmax] = bounds(squeeze(s_R{ii}(1,:,:)),'all');
    [ymin,ymax] = bounds(squeeze(s_R{ii}(2,:,:)),'all');
    [zmin,zmax] = bounds(squeeze(s_R{ii}(3,:,:)),'all');
    
    bx(ii,:) = round([floor(xmin-5),ceil(xmax+5)],-1);
    by(ii,:) = round([floor(ymin-5),ceil(ymax+5)],-1);
    bz(ii,:) = round([floor(zmin-5),ceil(zmax+5)],-1);
end

%% plot
fn = fn+1;
figure(fn)
set(gcf,'Position',[200 100 2*560 2*420])

% % % video setting
video = VideoWriter('vid_Test', 'Motion JPEG AVI');
video.FrameRate = 25*1/resampleDataRate;

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
        xlim([-max(abs(bx(:)))-10 max(abs(bx(:)))+10]);
        ylim([-max(abs(by(:))) max(abs(by(:)))]);
        zlim([0 max(bz(:))]);
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


