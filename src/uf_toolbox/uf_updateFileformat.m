function ufresult = uf_updateFileformat(ufresult)
% updates the ufresult file struct to some changes in variable names


ufresult = rename_stru(ufresult,'deconv','unfold'); 
ufresult.unfold = rename_stru(ufresult.unfold,'dcX','Xdc');
ufresult.unfold = rename_stru(ufresult.unfold,'dcX_termidx','Xdc_terms2cols');
ufresult.unfold = rename_stru(ufresult.unfold,'dcBeta','beta_dc');
ufresult.unfold = rename_stru(ufresult.unfold,'XBeta','beta_nodc');
ufresult.unfold = rename_stru(ufresult.unfold,'dcBasis','timebasis');
ufresult.unfold = rename_stru(ufresult.unfold,'dcX_termidx','Xdc_terms2cols');
ufresult.unfold = rename_stru(ufresult.unfold,'dcBasistime','times');
ufresult.unfold = rename_stru(ufresult.unfold,'col2eventtypes','cols2eventtypes');


ufresult.unfold = rename_stru(ufresult.unfold,'variableType','variabletypes');
ufresult.unfold = rename_stru(ufresult.unfold,'variableNames','variablenames');
ufresult.unfold = rename_stru(ufresult.unfold,'cols2variableNames','cols2variablenames');
ufresult.unfold = rename_stru(ufresult.unfold,'cols2eventtype','cols2eventtypes');
ufresult.unfold = rename_stru(ufresult.unfold,'eventtype','eventtypes');

ufresult.unfold = rename_stru(ufresult.unfold,'predictorSplines','splines');




ufresult = rename_stru(ufresult,'epoch','param');


end

function stru = rename_stru(stru,old,new)
if isfield(stru,old)
    stru.(new) = stru.(old);
    stru = rmfield(stru,old);
end
end
