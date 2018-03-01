function ufresult = uf_updateFileformat(ufresult)
% updates the ufresult file struct to some changes in variable names
unfold.deconv = rename_stru(ufresult.deconv,'dcX','Xdc');
unfold.deconv = rename_stru(ufresult.deconv,'dcX_termidx','Xuf_terms2cols');
unfold.deconv = rename_stru(ufresult.deconv,'dcBeta','beta_dc');
unfold.deconv = rename_stru(ufresult.deconv,'XBeta','beta_nodc');
unfold.deconv = rename_stru(ufresult.deconv,'dcBasis','timebasis');
unfold.deconv = rename_stru(ufresult.deconv,'dcX_termidx','Xuf_terms2cols');
unfold.deconv = rename_stru(ufresult.deconv,'dcBasistime','times');
unfold.deconv = rename_stru(ufresult.deconv,'col2eventtypes','cols2eventtypes');


unfold.deconv = rename_stru(ufresult.deconv,'variableType','variabletypes');
unfold.deconv = rename_stru(ufresult.deconv,'variableNames','variablenames');
unfold.deconv = rename_stru(ufresult.deconv,'cols2variableNames','cols2variablenames');
unfold.deconv = rename_stru(ufresult.deconv,'cols2eventtype','cols2eventtypes');
unfold.deconv = rename_stru(ufresult.deconv,'eventtype','eventtypes');

unfold.deconv = rename_stru(ufresult.deconv,'predictorSplines','splines');




ufresult = rename_stru(ufresult,'epoch','param');


end

function stru = rename_stru(stru,old,new)
if isfield(stru,old)
    stru.(new) = stru.(old);
    stru = rmfield(stru,old);
end
end