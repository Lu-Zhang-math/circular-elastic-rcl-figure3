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
cacheVersion=2;
residualTolerance=1e-10;
sourceHash=source_hash({[mfilename('fullpath'),'.m'], ...
    fullfile(root,'solve_elastic_rcl.m')});

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
                'sourceHash',sourceHash,'omega',omegaValues(iw), ...
                'L',LValues(iL),'params',params);
            useCache=false;
            if isfile(file)
                try
                    saved=load(file,'out','caseInfo');
                    useCache=isfield(saved,'out') && ...
                        isfield(saved,'caseInfo') && ...
                        isequaln(saved.caseInfo,caseInfo) && ...
                        valid_case_output(saved.out,residualTolerance);
                catch
                    useCache=false;
                end
            end
            if useCache
                out=saved.out;
            else
                out=solve_elastic_rcl(omegaValues(iw),LValues(iL),params);
                if ~valid_case_output(out,residualTolerance)
                    error('run_figure3:InvalidResult', ...
                        'Invalid result for d=%g, omega=%g, L=%g.', ...
                        dValues(id),omegaValues(iw),LValues(iL));
                end
                save_case(file,out,caseInfo);
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
fitIntercept=zeros(nd,nw);
fitR2=zeros(nd,nw);
for id=1:nd
    for iw=1:nw
        fitError=squeeze(errorH1(id,iw,fitMask)).';
        if any(~isfinite(fitError)) || any(fitError<=0)
            error('run_figure3:InvalidFitData', ...
                'Invalid fit data for d=%g, omega=%g.', ...
                dValues(id),omegaValues(iw));
        end
        fitL=LValues(fitMask);
        fitLogError=log(fitError);
        c=polyfit(fitL,fitLogError,1);
        fitted=polyval(c,fitL);
        total=sum((fitLogError-mean(fitLogError)).^2);
        residualFit=sum((fitLogError-fitted).^2);
        rate(id,iw)=-c(1);
        fitIntercept(id,iw)=c(2);
        fitR2(id,iw)=1-residualFit/max(total,realmin);
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
data=sortrows(data,{'d','omega','L'});
writetable(data,fullfile(output,'figure3_data.csv'));

[dRate,omegaRate]=ndgrid(dValues,omegaValues);
fitData=table(dRate(:),omegaRate(:),rate(:),fitIntercept(:),fitR2(:), ...
    'VariableNames',{'d','omega','fittedRate','intercept','R2'});
fitData=sortrows(fitData,{'d','omega'});
writetable(fitData,fullfile(output,'figure3_fitted_rates.csv'));

settings=struct('d',dValues,'omega',omegaValues,'L',LValues, ...
    'degree',degrees,'base',base,'cacheVersion',cacheVersion, ...
    'sourceHash',sourceHash,'fitRange',[4,9], ...
    'residualTolerance',residualTolerance, ...
    'matlabVersion',version,'matlabRelease',version('-release'));
results=struct('settings',settings,'errorH1',errorH1, ...
    'residual',residual,'backwardError',backwardError, ...
    'rate',rate,'fitIntercept',fitIntercept,'fitR2',fitR2, ...
    'data',data,'fitData',fitData);
save(fullfile(output,'figure3_results.mat'),'results');

fig=figure('Color','w','Visible','off','Position',[40,80,1500,450]);
layout=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
colors=lines(nw);
markers={'o','s','^'};
for id=1:nd
    ax=nexttile(layout,id);
    ax.YScale='log';
    hold(ax,'on');
    for iw=1:nw
        semilogy(ax,LValues,squeeze(errorH1(id,iw,:)), ...
            ['-',markers{iw}],'Color',colors(iw,:), ...
            'LineWidth',1.2,'MarkerSize',4.5,'MarkerFaceColor','w', ...
            'DisplayName',sprintf('$\\omega=%g\\,(p=%.2f)$', ...
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

function tf=valid_case_output(out,residualTolerance)
required={'omega','L','uRelH1','scaledResidual','backwardError'};
tf=isstruct(out) && all(isfield(out,required));
if ~tf
    return;
end
values=[out.omega,out.L,out.uRelH1,out.scaledResidual,out.backwardError];
tf=all(isfinite(values)) && out.uRelH1>0 && ...
    out.scaledResidual>=0 && out.backwardError>=0 && ...
    out.scaledResidual<=residualTolerance && ...
    out.backwardError<=residualTolerance;
end

function save_case(file,out,caseInfo)
folder=fileparts(file);
temporary=[tempname(folder),'.mat'];
cleanup=onCleanup(@() delete_if_present(temporary));
save(temporary,'out','caseInfo');
[ok,message]=movefile(temporary,file,'f');
if ~ok
    error('run_figure3:CacheWriteFailed','%s',message);
end
clear cleanup
end

function delete_if_present(file)
if isfile(file)
    delete(file);
end
end

function value=source_hash(files)
digest=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
for k=1:numel(files)
    fid=fopen(files{k},'rb');
    if fid<0
        error('run_figure3:SourceReadFailed', ...
            'Cannot read source file: %s',files{k});
    end
    cleanup=onCleanup(@() fclose(fid));
    bytes=fread(fid,Inf,'*uint8');
    digest.update(typecast([bytes;uint8(0)],'int8'));
    clear cleanup
end
raw=typecast(digest.digest(),'uint8');
value=lower(reshape(dec2hex(raw,2).',1,[]));
end
