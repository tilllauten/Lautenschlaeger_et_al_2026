import colorcet as cc

"""Collection of colours from various sources aiming to be colour-blind friendly."""

# https://www.nature.com/articles/nmeth.1618
wong_colours = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

# https://personal.sron.nl/~pault/
tol_bright = ['#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']
tol_vibrant = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB']
tol_muted = ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499']
tol_highcontrast = ['#004488', '#DDAA33', '#BB5566']

# https://davidmathlogic.com/colorblind/
two_contrasts = [['#E1BE6A', '#40B0A6'],  # orange, teal
                 ['#FFC20A', '#0C7BDC'],  # yellow blue
                 ['#994F00', '#006CD1']]  # brown blue

# Longer list of colours for categorical colours https://colorcet.holoviz.org/user_guide/Categorical.html
# Generated based on Glasbey, Chris; van der Heijden, Gerie & Toh, Vivian F. K. et al. (2007),
# “Colour displays for categorical images”, Color Research & Application 32.4: 304-309.
categ_colours = cc.glasbey_bw_minc_20
del categ_colours[19]  # The yellow colour is almost invisible.
glasbey_palettes = {'glasbey': categ_colours,
                    'glasbey_cool': cc.glasbey_cool,
                    'glasbey_warm': cc.glasbey_warm,
                    'glasbey_dark': cc.glasbey_dark,
                    'glasbey_light': cc.glasbey_light}


# Own palettes
red_colours = ['#2e0000', '#a11b1b', '#ed6109', '#ffd86b', '#ffc400']
avghic = ['#2e0000', '#ab153f', '#1231c9', '#61e8e8']
fulco_colours = ['#387329', '#afcf65', '#1231c9', '#46cdf2', '#9bfa02']
poster_colours = ['#00003b', '#1231c9', '#46cdf2', '#d6c209']
vs_hic = ['#a11b1b', "#f2833d", '#033e61', '#84d7e3'] * 2
enformer = ['#a11b1b', '#033e61', '#bbcc33', "#f5f069", '#a11b1b', "#f2833d"] * 2
blue_colours = ['#033e61', '#6445ff', '#08d19c', '#95ff8a', '#ffc400']
gabc_avghic = ['#033e61', '#84d7e3', '#fab30f'] * 2

# CARE the following shapes are not plotted for whatever reason: ['1', '2', '3', '4', '+', 'x']
marker_shapes = ['o', 'D', 'P', 's', '*', 'X', '<', 'p', 'd', '^', 'v', '>', 'H', '$O$', '$D$', '$U$', '$Y$', '$N$']

hatches = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']


