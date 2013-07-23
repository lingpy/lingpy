from ...settings import rcParams

phybo = dict(
                phybo_figsize          = (10,10),
                phybo_linescale        = 1.0,
                phybo_maxweight        = False,
                phybo_xlim             = 5,
                phybo_ylim             = 5,
                phybo_xlimr            = False,
                phybo_xliml            = False,
                phybo_ylimt            = False,
                phybo_ylimb            = False,
                phybo_left             = 0.01,
                phybo_right            = 0.99,
                phybo_top              = 0.99,
                phybo_bottom           = 0.01,
                phybo_cbar_shrink      = 0.55,
                phybo_cbar_fraction    = 0.1,
                phybo_cbar_pad         = 0.1,
                phybo_cbar_orientation = 'vertical',
                phybo_cbar_label       = 'Inferred Links',
                phybo_vedgestyle       = 'double',
                phybo_vedgecolor       = 'black',
                phybo_vedgelinewidth   = 5,
                phybo_hedgescale       = 3,
                phybo_nodestyle        = 'double',
                phybo_nodesize         = 10,
                phybo_nodecolor        = 'black',
                phybo_labels           = {},
                phybo_prefix = '- ',
                phybo_suffix = ' -',
                phybo_textsize = '10',
                phybo_vsd_scale = 0.1,
                phybo_latex_preamble = [],
                phybo_fileformat = "pdf"
                )
rcParams.update(phybo)