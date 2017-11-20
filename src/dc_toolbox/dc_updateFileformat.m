function unfold = dc_updateFileformat(unfold)
% updates the unfold file struct to some changes in variable names
unfold.deconv = rename_deconv(unfold.deconv,'dcX','Xdc');
unfold.deconv = rename_deconv(unfold.deconv,'dcX_termidx','Xdc_terms2cols');
unfold.deconv = rename_deconv(unfold.deconv,'dcBeta','beta_dc');
unfold.deconv = rename_deconv(unfold.deconv,'XBeta','beta_nodc');
unfold.deconv = rename_deconv(unfold.deconv,'dcBasis','timebasis');
unfold.deconv = rename_deconv(unfold.deconv,'dcX_termidx','Xdc_terms2cols');
unfold.deconv = rename_deconv(unfold.deconv,'dcBasistime','times');

end
function deconv = rename_deconv(deconv,old,new)
if isfield(deconv,old)
    deconv.(new) = deconv.(old);
    deconv = rmfield(deconv,old);
end
end