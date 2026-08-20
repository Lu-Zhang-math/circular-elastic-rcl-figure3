function results=run_figure3

root=fileparts(mfilename('fullpath'));
output=fullfile(root,'results');
cases=fullfile(output,'cases');
if ~isfolder(cases)
    mkdir(cases);
end

dValues=[0.01,0.1,1];
omegaValues=[5,10,50];
LValues=2:14;
degrees=[14,40,24;18,45,30;30,60,40];
base=struct('rho',2,'lambda',3,'mu',1,'r0',0.5,'a',1, ...
    'source',[0.2,0.1],'M',10,'nPhysicalLayers',4);
cacheVersion=1;

nd=numel(dValues);
nw=numel(omegaValues);
nL=numel(LValues);
errorH1=zeros(nd,nw,nL);
residual=zeros(nd,nw,nL);
backwardError=zeros(nd,nw,nL);

for id=1:nd
    for iw=1:nw
        params=base;
        params.b=params.a+dValues(id);
        params.N1=degrees(iw,1);
        params.N2=degrees(iw,2);
        params.N3=degrees(iw,3);
        for iL=1:nL
            key=sprintf('d_%s_w_%g_L_%g', ...
                strrep(sprintf('%g',dValues(id)),'.','p'), ...
                omegaValues(iw),LValues(iL));
            file=fullfile(cases,[key,'.mat']);
            caseInfo=struct('version',cacheVersion, ...
                'omega',omegaValues(iw),'L',LValues(iL),'params',params);
            useCache=false;
            if isfile(file)
                saved=load(file,'out','caseInfo');
                useCache=isfield(saved,'out') && ...
                    isfield(saved,'caseInfo') && ...
                    isequaln(saved.caseInfo,caseInfo);
            end
            if useCache
                out=saved.out;
            else
                out=solve_elastic_rcl(omegaValues(iw),LValues(iL),params);
                save(file,'out','caseInfo');
            end
            errorH1(id,iw,iL)=out.uRelH1;
            residual(id,iw,iL)=out.scaledResidual;
            backwardError(id,iw,iL)=out.backwardError;
            fprintf('d=%g, omega=%g, L=%g: %.6e\n', ...
                dValues(id),omegaValues(iw),LValues(iL),out.uRelH1);
        end
    end
end

fitMask=LValues>=4 & LValues<=9;
rate=zeros(nd,nw);
for id=1:nd
    for iw=1:nw
        c=polyfit(LValues(fitMask), ...
            log(squeeze(errorH1(id,iw,fitMask))).',1);
        rate(id,iw)=-c(1);
    end
end

[dGrid,omegaGrid,LGrid]=ndgrid(dValues,omegaValues,LValues);
N1Grid=repmat(reshape(degrees(:,1),1,nw,1),nd,1,nL);
N2Grid=repmat(reshape(degrees(:,2),1,nw,1),nd,1,nL);
N3Grid=repmat(reshape(degrees(:,3),1,nw,1),nd,1,nL);
data=table(dGrid(:),omegaGrid(:),LGrid(:),N1Grid(:),N2Grid(:),N3Grid(:), ...
    errorH1(:),residual(:),backwardError(:), ...
    'VariableNames',{'d','omega','L','N1','N2','N3', ...
    'relativeH1Error','scaledResidual','backwardError'});
writetable(data,fullfile(output,'figure3_data.csv'));

[dRate,omegaRate]=ndgrid(dValues,omegaValues);
fitData=table(dRate(:),omegaRate(:),rate(:), ...
    'VariableNames',{'d','omega','fittedRate'});
writetable(fitData,fullfile(output,'figure3_fitted_rates.csv'));

settings=struct('d',dValues,'omega',omegaValues,'L',LValues, ...
    'degree',degrees,'base',base,'cacheVersion',cacheVersion);
results=struct('settings',settings,'errorH1',errorH1, ...
    'rate',rate,'data',data);
save(fullfile(output,'figure3_results.mat'),'results');

fig=figure('Color','w','Visible','off','Position',[40,80,1500,450]);
layout=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
colors=lines(nw);
markers={'o','s','^'};
for id=1:nd
    ax=nexttile(layout,id);
    hold(ax,'on');
    for iw=1:nw
        semilogy(ax,LValues,squeeze(errorH1(id,iw,:)), ...
            ['-',markers{iw}],'Color',colors(iw,:), ...
            'LineWidth',1.2,'MarkerSize',4.5,'MarkerFaceColor','w', ...
            'DisplayName',sprintf('$\\omega=%g\,(p=%.2f)$', ...
            omegaValues(iw),rate(id,iw)));
    end
    reference=errorH1(id,1,1)*exp(-(LValues-LValues(1)));
    semilogy(ax,LValues,reference,'k--','LineWidth',1.2, ...
        'DisplayName','$C\mathrm{e}^{-L}$');
    grid(ax,'on');
    box(ax,'on');
    xlim(ax,[LValues(1),LValues(end)]);
    xticks(ax,2:2:14);
    xlabel(ax,'$L=\tau_0d$','Interpreter','latex');
    ylabel(ax,'relative $H^1(\Omega_1)$ error','Interpreter','latex');
    title(ax,sprintf('$d=%g$',dValues(id)),'Interpreter','latex');
    legend(ax,'Location','southwest','Interpreter','latex');
    ax.TickLabelInterpreter='latex';
    ax.FontName='Times New Roman';
    ax.FontSize=11;
end
exportgraphics(fig,fullfile(output,'figure3.png'),'Resolution',300);
exportgraphics(fig,fullfile(output,'figure3.pdf'),'ContentType','vector');
savefig(fig,fullfile(output,'figure3.fig'));
close(fig);
end
