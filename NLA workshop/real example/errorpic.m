%Generates the error figure
for j=1:2
plot(log(H),log(errmatrix(:,j)),'o-');
hold on
end