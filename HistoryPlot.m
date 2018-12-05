DV = zeros(6400,80);
for ii = 1:1:80
    FileName = ['.\DesignVariables\','DV_Iter',num2str(ii),'.dat'];
    DV(:,ii) = load(FileName);
end