"""
Specific settings for PhyBo class.
"""
from __future__ import unicode_literals, division, print_function

from lingpy.settings import rcParams

rcParams.update(
    phybo_bottom=0.01,
    phybo_cbar_fraction=0.1,
    phybo_cbar_label='Inferred Links',
    phybo_cbar_orientation='vertical',
    phybo_cbar_pad=0.1,
    phybo_cbar_shrink=0.55,
    phybo_figsize=(10, 10),
    phybo_fileformat="pdf",
    phybo_hedgescale=3,
    phybo_labels={},
    phybo_latex_preamble=[],
    phybo_left=0.01,
    phybo_linescale=1.0,
    phybo_maxweight=False,
    phybo_nodecolor='black',
    phybo_nodesize=10,
    phybo_nodestyle='double',
    phybo_prefix='- ',
    phybo_right=0.99,
    phybo_suffix=' -',
    phybo_textsize='10',
    phybo_top=0.99,
    phybo_vedgecolor='black',
    phybo_vedgelinewidth=5,
    phybo_vedgeinnerline=3,
    phybo_vedgestyle='double',
    phybo_vsd_scale=0.1,
    phybo_xlim=5,
    phybo_xliml=False,
    phybo_xlimr=False,
    phybo_ylim=5,
    phybo_ylimb=False,
    phybo_ylimt=False,
)
