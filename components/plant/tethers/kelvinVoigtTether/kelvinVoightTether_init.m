vssBlk = [gcb, '/kelvinVoigtTether'];

oldVariant = get_param(vssBlk,'OverrideUsingVariant');

if any([class_thr.numNodes] == 2)
    set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNode2_variant');
else
    set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNodeN_variant');
end