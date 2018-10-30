function Wdiff = numdiff(W, direct, dx)
Wdiff = zeros(size(W));
if direct == 1
    Wdiff(2:end-1,:) = (W(3:end,:) - W(1:end-2,:))/(2*dx);
    Wdiff(1,:) = 2*Wdiff(2,:) - Wdiff(3,:);
    Wdiff(end,:) = 2*Wdiff(end-1,:) - Wdiff(end-2,:);
elseif direct == 2
    Wdiff(:,2:end-1) = (W(:,3:end) - W(:,1:end-2))/(2*dx);
    Wdiff(:,1) = 2*Wdiff(:,2) - Wdiff(:,3);
    Wdiff(:,end) = 2*Wdiff(:,end-1) - Wdiff(:,end-2);
else
    disp('wrong direction')
end
end

% Test
% Test
% whaddap
testnumdiff = 400;