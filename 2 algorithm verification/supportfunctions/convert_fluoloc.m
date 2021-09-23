function [ x, xXYZind, xAind ] = convert_fluoloc(X,Y,Z,XYZloc)
%Convert XYZloc in um to logical indices in X,Y,Z

XYZloc(XYZloc(:,1) < X(1),1) =  X(1);
XYZloc(XYZloc(:,1) > X(end),1) =  X(end);
XYZloc(XYZloc(:,2) < Y(1),2) =  Y(1);
XYZloc(XYZloc(:,2) > Y(end),2) =  Y(end);
XYZloc(XYZloc(:,3) < Z(1),3) =  Z(1);
XYZloc(XYZloc(:,3) > Z(end),3) =  Z(end);


Xsize = length(X);
Ysize = length(Y);
Zsize = length(Z);
if length(Z) == 1
    stepsize = [X(2)-X(1) Y(2)-Y(1) Z(end)-Z(1)];
else
    stepsize = [X(2)-X(1) Y(2)-Y(1) Z(2)-Z(1)];
end
x = zeros(Xsize*Ysize*Zsize,1); %this is the indexed x-vector y=Ax

XYZlocfromorig = bsxfun(@minus,XYZloc,[X(1) Y(1) Z(1)]);
xXYZind = round(bsxfun(@rdivide,XYZlocfromorig,stepsize)+1); %require the offset

xXYZind(isnan(xXYZind)) = 1;
xAind = sub2ind([Xsize Ysize Zsize],xXYZind(:,1),xXYZind(:,2),xXYZind(:,3));

% [~,I] = sort(xAind,'ascend');
% xAind = xAind(I);
% xXYZind = xXYZind(I,:);
x(xAind) = 1;

end