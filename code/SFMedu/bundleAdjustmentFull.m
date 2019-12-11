function graph = bundleAdjustmentFull(graph)

% convert from Rt matrix to AngleAxis
nCam=length(graph.frames);
Mot = zeros(3,2,nCam);
for camera=1:nCam
    Mot(:,1,camera) = RotationMatrix2AngleAxis(graph.Mot(:,1:3,camera));
    Mot(:,2,camera) = graph.Mot(:,4,camera);
end



Str = graph.Str;
f  = graph.f;

if ~isfield(graph, 'fx') || ~isfield(graph, 'fy') || ~isfield(graph, 'px') || ~isfield(graph, 'py')
    fx = f;
    fy = f;
    px = 0;
    py = 0;
else
    fx = graph.fx;
    fy = graph.fy;
    px = graph.px;
    py = graph.py;
end

residuals = reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,Str);
fprintf('initial error = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));

% bundle adjustment using lsqnonlin in Matlab (Levenberg-Marquardt)
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','off');

% adjust structure [for homework]
% !!! fill in your code here
[vec,resnorm,residuals,exitflag] = lsqnonlin(@(optme) reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,optme, Str), [Mot(:)],[],[],options);

% reshape optimized result 
Mot = reshape(vec(1:6*nCam), 3, 2, []);
fprintf('motion only error = %f\n', 2*sqrt(resnorm/length(residuals)));

% adjust motion [for homework]
[vec,resnorm,residuals,exitflag] = lsqnonlin(@(optme) reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,optme), [Str(:)],[],[],options);
% reshape optimized result\
Str = reshape(vec, 3, []);
% !!! fill in your code here
fprintf('structure only error = %f\n', 2*sqrt(resnorm/length(residuals)));

% adjust motion and structure
[vec,resnorm,residuals,exitflag] = lsqnonlin(@(x) reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,x), [Mot(:); Str(:)],[],[],options);
% [Mot,Str] = unpackMotStrf(nCam,vec);
fprintf('S&M error = %f\n', 2*sqrt(resnorm/length(residuals)));

[vec,resnorm,residuals,exitflag] = lsqnonlin(@(x) reprojectionResidualFull(graph.ObsIdx,graph.ObsVal,px,py,fx,fy,Mot,Str), [fx; fy; px; py; Mot(:); Str(:)],[],[],options);
[Mot,Str,fxnew,fynew,pxnew,pynew] = unpackMotStrFull(nCam,vec);
fprintf('S&M error with adjust focal length = %f\n', resnorm/length(residuals));
graph.f = fxnew;
graph.fx = fxnew;
graph.fy = fynew;
graph.px = pxnew;
graph.py = pynew;

%residuals = reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,Str);
%fprintf('final error = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));


for camera=1:nCam
    graph.Mot(:,:,camera) = [AngleAxis2RotationMatrix(Mot(:,1,camera))  Mot(:,2,camera)];    
end
graph.Str = Str;

