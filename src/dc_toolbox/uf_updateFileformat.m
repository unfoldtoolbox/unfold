function ufresult = uf_updateFileformat(ufresult)
% updates the ufresult file struct to some changes in variable names
ufresult = rename_stru(ufresult,'deconv','unfold'); 
ufresult.deconv = rename_stru(ufresult.deconv,'dcX','Xdc');
ufresult.deconv = rename_stru(ufresult.deconv,'dcX_termidx','Xuf_terms2cols');
ufresult.deconv = rename_stru(ufresult.deconv,'dcBeta','beta_dc');
ufresult.deconv = rename_stru(ufresult.deconv,'XBeta','beta_nodc');
ufresult.deconv = rename_stru(ufresult.deconv,'dcBasis','timebasis');
ufresult.deconv = rename_stru(ufresult.deconv,'dcX_termidx','Xuf_terms2cols');
ufresult.deconv = rename_stru(ufresult.deconv,'dcBasistime','times');
ufresult.deconv = rename_stru(ufresult.deconv,'col2eventtypes','cols2eventtypes');


ufresult.deconv = rename_stru(ufresult.deconv,'variableType','variabletypes');
ufresult.deconv = rename_stru(ufresult.deconv,'variableNames','variablenames');
ufresult.deconv = rename_stru(ufresult.deconv,'cols2variableNames','cols2variablenames');
ufresult.deconv = rename_stru(ufresult.deconv,'cols2eventtype','cols2eventtypes');
ufresult.deconv = rename_stru(ufresult.deconv,'eventtype','eventtypes');

ufresult.deconv = rename_stru(ufresult.deconv,'predictorSplines','splines');




ufresult = rename_stru(ufresult,'epoch','param');


end

function stru = rename_stru(stru,old,new)
if isfield(stru,old)
    stru.(new) = stru.(old);
    stru = rmfield(stru,old);
end
end
