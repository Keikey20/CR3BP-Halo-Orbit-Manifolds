
outdir = 'halo_output';
load(fullfile(outdir,'monodromy.mat'),'eigvals','eigvecs','PhiT'); %#ok<LOAD>
load(fullfile(outdir,'X0_corrected.mat'),'X0','T','mu'); %#ok<LOAD>
load(fullfile(outdir,'manifolds.mat'),'manif'); %#ok<LOAD>
%compute collinear Lagrange points% 
options = optimset('Display','off');
f = @(x) dOmega_dx_on_x(x, mu);
% initial guesses L1 & L2
xL1 = fsolve(f, 0.99, options);
xL2 = fsolve(f, 1.02, options);
% Annotated 3d%
fig = figure('Name','Annotated 3D manifold','Position',[50 50 900 700]); hold on; grid on; axis equal
plot3(-mu,0,0,'ko','MarkerFaceColor','k','MarkerSize',8);
text(-mu,0,0,'  Primary 1','FontSize',10);
plot3(1-mu,0,0,'ro','MarkerFaceColor','r','MarkerSize',8);
text(1-mu,0,0,'  Primary 2','FontSize',10);
% plot halo orb
[tpo, Ypo] = ode45(@(t,y) cr3bp_EOM(t,y,mu), linspace(0,T,1000), X0);
plot3(Ypo(:,1),Ypo(:,2),Ypo(:,3),'k-','LineWidth',1.8,'DisplayName','Halo');
% markin' L1 and L2
plot3(xL1,0,0,'mp','MarkerSize',12,'MarkerFaceColor','m');
text(xL1,0,0,' L1','Color','m','FontSize',10);
plot3(xL2,0,0,'md','MarkerSize',12,'MarkerFaceColor','m');
text(xL2,0,0,' L2','Color','m','FontSize',10);
% plot manifolds with dif colors
for k=1:length(manif)
    Y = manif(k).Y;
    if strcmp(manif(k).type,'unstable')
        plot3(Y(:,1),Y(:,2),Y(:,3),'-','LineWidth',1.2,'Color',[0.85 0.3 0.3]);
    else
        plot3(Y(:,1),Y(:,2),Y(:,3),'--','LineWidth',1.2,'Color',[0.2 0.4 0.8]);
    end
    plot3(manif(k).X0(1),manif(k).X0(2),manif(k).X0(3),'o','MarkerFaceColor',[0 0 0],'MarkerSize',4);
end

xlabel('x'); ylabel('y'); zlabel('z'); title('Halo orbit and annotated invariant manifolds');
legend({'Primary1','Primary2','Halo','L-points'}, 'Location','bestoutside');
view(3);
saveas(fig, fullfile(outdir,'annotated_manifolds_3D.png'));
% XY, XZ, YZ projections
projNames = {'xy','xz','yz'}; idxs = {[1 2],[1 3],[2 3]};
for p=1:3
    fh = figure('Name',['Projection ' projNames{p}],'Position',[50 50 700 600]); hold on; grid on
    plot(Ypo(:,idxs{p}(1)), Ypo(:,idxs{p}(2)),'k-','LineWidth',1.6);
    plot(-mu,0,'ko','MarkerFaceColor','k'); plot(1-mu,0,'ro','MarkerFaceColor','r');
    for k=1:length(manif)
        Y = manif(k).Y;
        if strcmp(manif(k).type,'unstable')
            plot(Y(:,idxs{p}(1)), Y(:,idxs{p}(2)),'-','LineWidth',1.1,'Color',[0.85 0.3 0.3]);
        else
            plot(Y(:,idxs{p}(1)), Y(:,idxs{p}(2)),'--','LineWidth',1.1,'Color',[0.2 0.4 0.8]);
        end
    end
    xlabel(projNames{p}(1)); ylabel(projNames{p}(2));
    title(['Projection: ' upper(projNames{p})]);
    saveas(fh, fullfile(outdir,['projection_' projNames{p} '.png']));
end
% Refining Poincare crossings%
pc_un = []; pc_st = [];
for k=1:length(manif)
    Y = manif(k).Y;
    tt = manif(k).t;
   
    s = sign(Y(:,2));
    idx = find(s(1:end-1).*s(2:end) < 0);
    for ii = idx'
        % linear interpolation between ii and ii+1
        alpha = abs(Y(ii,2)) / (abs(Y(ii,2))+abs(Y(ii+1,2)));
        xcross = Y(ii,1) + alpha*(Y(ii+1,1)-Y(ii,1));
        vxcross = Y(ii,4) + alpha*(Y(ii+1,4)-Y(ii,4));
        if strcmp(manif(k).type,'unstable')
            pc_un(end+1,:) = [xcross, vxcross]; %#ok<SAGROW>
        else
            pc_st(end+1,:) = [xcross, vxcross]; %#ok<SAGROW>
        end
    end
end
% plot Poincare map
figp = figure('Name','Poincare Map (y=0)','Position',[50 50 700 600]); hold on; grid on;
if ~isempty(pc_un), scatter(pc_un(:,1), pc_un(:,2), 40, [0.85 0.3 0.3], 'filled', 'DisplayName','Unstable crossings'); end
if ~isempty(pc_st), scatter(pc_st(:,1), pc_st(:,2), 40, [0.2 0.4 0.8], 'filled', 'DisplayName','Stable crossings'); end
xlabel('x'); ylabel('vx'); title('Poincaré section (y=0): x vs v_x'); legend('Location','best');
saveas(figp, fullfile(outdir,'poincare_colored.png'));

fprintf('Annotated figures saved in folder: %s\n', outdir);
%% derivative dOmega/dx on x-axis y=z=0
function val = dOmega_dx_on_x(x, mu)
    r1 = abs(x+mu);
    r2 = abs(x-1+mu);
    val = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
end
