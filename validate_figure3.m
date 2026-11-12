function report=validate_figure3(results)

if nargin==0
    root=fileparts(mfilename('fullpath'));
    file=fullfile(root,'results','figure3_results.mat');
    loaded=load(file,'results');
    results=loaded.results;
end

expectedSize=[3,3,13];
if ~isequal(size(results.errorH1),expectedSize) || ...
        ~isequal(size(results.residual),expectedSize) || ...
        ~isequal(size(results.backwardError),expectedSize)
    error('validate_figure3:InvalidSize','Unexpected result dimensions.');
end

if any(~isfinite(results.errorH1),'all') || ...
        any(results.errorH1<=0,'all') || ...
        any(~isfinite(results.residual),'all') || ...
        any(~isfinite(results.backwardError),'all')
    error('validate_figure3:InvalidData','Results contain invalid values.');
end

residualTolerance=results.settings.residualTolerance;
maxResidual=max(results.residual,[],'all');
maxBackwardError=max(results.backwardError,[],'all');
if maxResidual>residualTolerance || maxBackwardError>residualTolerance
    error('validate_figure3:ResidualTooLarge', ...
        'Residual tolerance was not met.');
end

expectedRates=[1.54,1.53,1.51;1.54,1.53,1.51;1.53,1.52,1.51];
rateTolerance=0.03;
maxRateDifference=max(abs(results.rate-expectedRates),[],'all');
if maxRateDifference>rateTolerance
    error('validate_figure3:RateMismatch', ...
        'Fitted rates do not match the reference values.');
end

id=find(abs(results.settings.d-0.1)<1e-12,1);
iw=find(abs(results.settings.omega-5)<1e-12,1);
iL=find(abs(results.settings.L-2)<1e-12,1);
referenceError=6.3015954052119856e-3;
representativeError=results.errorH1(id,iw,iL);
relativeDifference=abs(representativeError-referenceError)/referenceError;
if relativeDifference>1e-4
    error('validate_figure3:ReferenceMismatch', ...
        'The representative error does not match the reference value.');
end

report=struct('passed',true,'maxResidual',maxResidual, ...
    'maxBackwardError',maxBackwardError, ...
    'maxRateDifference',maxRateDifference, ...
    'representativeError',representativeError);

fprintf('Figure 3 validation passed.\n');
fprintf('Maximum scaled residual: %.3e\n',maxResidual);
fprintf('Maximum backward error: %.3e\n',maxBackwardError);
fprintf('Maximum fitted-rate difference: %.3e\n',maxRateDifference);
end
