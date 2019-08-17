function [objTbl,obj] = objectiveFunction(X)

x1 = X(:,1);
x2 = X(:,2);

obj = 1*(-((x1.^2 + x2.^2)./50) + 1);

objTbl = table(obj(:),x1(:),x2(:),'VariableNames',{'objF','x1','x2'});

end