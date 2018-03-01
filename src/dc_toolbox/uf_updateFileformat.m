function unfold = dc_updateFileformat(unfold)
% updates the unfold file struct to some changes in variable names
unfold.deconv = rename_stru(unfold.deconv,'dcX','Xdc');
unfold.deconv = rename_stru(unfold.deconv,'dcX_termidx','Xdc_terms2cols');
unfold.deconv = rename_stru(unfold.deconv,'dcBeta','beta_dc');
unfold.deconv = rename_stru(unfold.deconv,'XBeta','beta_nodc');
unfold.deconv = rename_stru(unfold.deconv,'dcBasis','timebasis');
unfold.deconv = rename_stru(unfold.deconv,'dcX_termidx','Xdc_terms2cols');
unfold.deconv = rename_stru(unfold.deconv,'dcBasistime','times');
unfold.deconv = rename_stru(unfold.deconv,'col2eventtypes','cols2eventtypes');


unfold.deconv = rename_stru(unfold.deconv,'variableType','variabletypes');
unfold.deconv = rename_stru(unfold.deconv,'variableNames','variablenames');
unfold.deconv = rename_stru(unfold.deconv,'cols2variableNames','cols2variablenames');
unfold.deconv = rename_stru(unfold.deconv,'cols2eventtype','cols2eventtypes');
unfold.deconv = rename_stru(unfold.deconv,'eventtype','eventtypes');

unfold.deconv = rename_stru(unfold.deconv,'predictorSplines','splines');




unfold = rename_stru(unfold,'epoch','param');


end

function stru = rename_stru(stru,old,new)
if isfield(stru,old)
    stru.(new) = stru.(old);
    stru = rmfield(stru,old);
end
end