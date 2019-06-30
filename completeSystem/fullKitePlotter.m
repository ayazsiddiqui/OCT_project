dsgnData = tp.aeroDesignData;


figure

sol_outline = cell(5,1);

int_mat1 = zeros(5,3);
int_mat2 = zeros(5,3,n_steps);

for jj = 1:5
    
    for ii = 1:n_steps
        
        for kk = 1:5
            int_mat1(kk,:) = ( sol_Rcm_o(:,ii) + ...
                rotation_sequence(sol_euler(:,ii))*...
                cross(sol_OwB(:,ii),dsgnData.outlines(jj).pts(kk,:))' )';
        end
        int_mat2(:,:,ii) = int_mat1;
    end
    sol_outline{jj} = int_mat2;
    
end

figure(1)
for ii = 1:n_steps
    
    if ii > 1
        h = findall(gca,'type','line','color',red,'-or','color',black);
        delete(h);
    end
    
    for jj = 1:5
        p_kite = plot3(sol_outline{jj}(:,1,ii),...
        sol_outline{jj}(:,2,ii),...
        sol_outline{jj}(:,3,ii),...
        '-','linewidth',line_wd,'color',black);
    hold on
    grid on
    
    end
    
    waitforbuttonpress
    
    
end


