from collections import namedtuple
from dataclasses import dataclass


# 1. Define the structure for a single color pair
ColorPair = namedtuple('ColorPair', ['base', 'alternative'])

# 2. Define the Palette to hold them all
@dataclass
class Palette:
    red: ColorPair
    blue: ColorPair
    green: ColorPair
    purple: ColorPair
    orange: ColorPair
    cyan:   ColorPair
    berry: ColorPair
    gold: ColorPair

    black: ColorPair
    brown: ColorPair

# 3. Create the instance with your "Easy on the Eyes" colors
colors = Palette(
    red   = ColorPair(base='tab:red',   alternative='maroon'),
    blue  = ColorPair(base='tab:blue',  alternative='midnightblue'),
    green = ColorPair(base='tab:green', alternative='darkgreen'),
    purple = ColorPair(base='tab:purple', alternative='rebeccapurple'),
    orange = ColorPair(base='tab:orange', alternative='saddlebrown'),
    cyan   = ColorPair(base='tab:cyan',   alternative='teal'),
    berry = ColorPair(base='deeppink',   alternative='darkmagenta'),
    gold = ColorPair(base='goldenrod',   alternative='olive'),

    black = ColorPair(base='black',     alternative='dimgray'),
    brown = ColorPair(base='tab:brown', alternative='sienna')
)