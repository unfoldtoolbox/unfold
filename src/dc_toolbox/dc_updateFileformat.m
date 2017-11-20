function unfold = dc_updateFileformat(unfold)
% updates the unfold file struct to some changes in variable names
unfold.deconv = rename_stru(unfold.deconv,'dcX','Xdc');
unfold.deconv = rename_stru(unfold.deconv,'dcX_termidx','Xdc_terms2cols');
unfold.deconv = rename_stru(unfold.deconv,'dcBeta','beta_dc');
unfold.deconv = rename_stru(unfold.deconv,'XBeta','beta_nodc');
unfold.deconv = rename_stru(unfold.deconv,'dcBasis','timebasis');
unfold.deconv = rename_stru(unfold.deconv,'dcX_termidx','Xdc_terms2cols');
unfold.deconv = rename_stru(unfold.deconv,'dcBasistime','times');
unfold.deconv = rename_stru(unfold.deconv,'col2eventtype','cols2eventtype');



unfold = rename_stru(unfold,'epoch','param');


end

function stru = rename_stru(stru,old,new)
if isfield(stru,old)
    stru.(new) = stru.(old);
    stru = rmfield(stru,old);
end
end