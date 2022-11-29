#!/usr/bin/env python

from colour_aap import aa_colors, normalize_color, normalize_aesnn3_color
import cairo
import sys


WIDTH = 18
HEIGHT = 14
OFF_X = 2
OFF_Y = 18

normalize_aesnn3_color()

surface = cairo.PSSurface( sys.argv[1], WIDTH * 21 + 2 * OFF_X, HEIGHT + 2.4 * OFF_Y )

c = cairo.Context(surface)

if True:
    # clear background
    c.set_antialias( cairo.ANTIALIAS_NONE )
    c.set_line_cap(cairo.LINE_CAP_BUTT)
    c.set_source_rgb( 1, 1, 1 )
    c.set_operator(cairo.OPERATOR_SOURCE)
    c.paint()

off_x = OFF_X
off_y = OFF_Y
c.set_line_width(0.1)
c.set_font_size( HEIGHT - 4 )
fascent, fdescent, fheight, fxadvance, fyadvance = c.font_extents()
off_y_font = off_y + (HEIGHT + fheight)/2 - 1.25

for aa in 'CILVMAGFYWTHKRQNSDEPX':

    c.rectangle(off_x, off_y, WIDTH, HEIGHT)
    c.set_source_rgb( *aa_colors[aa] )
    c.fill()
    c.stroke()
    c.set_source_rgb( 0,0,0 )
    c.rectangle(off_x, off_y, WIDTH, HEIGHT)
    c.stroke()
    color_sum = sum( aa_colors[aa] )
    if color_sum > 1.25:
        c.set_source_rgb( 0, 0, 0)
    else:
        c.set_source_rgb( 1, 1, 1)
    xbearing, ybearing, width, height, xadv, yadv = (c.text_extents(aa))
    c.move_to( off_x + (WIDTH - xbearing - width)/2, off_y_font )
    c.show_text( aa )
    c.stroke()

    off_x += WIDTH

# top bar

upper_bars = [[ ('aliphatic', 2, 4), ('tiny', 6, 7), ('aromatic', 8, 10) ], [ ('hydrophobic', 1, 13), ('hydrophilic', 14, 20) ] ]
lower_bars = [ [ ('(+)', 12, 14), ('(-)', 18, 19)], [ ('non-polar', 1, 8), ('polar', 9, 19)] ] 

def draw_annotation( bars, offset, c, lower=True):
    c.set_source_rgb( 0,0,0 )
    c.set_line_width( 0.1 )
    c.set_font_size( 8 )
    fascent, fdescent, fheight, fxadvance, fyadvance = c.font_extents()
    for bar in bars:
        for (label, start_box, end_box) in bar:
            start_x = OFF_X + (start_box-1) * WIDTH + 2
            end_x = OFF_X + (end_box) * WIDTH - 2
            line_length = end_x - start_x
            xbearing, ybearing, width, height, xadv, yadv = (c.text_extents(label))
            if not lower:
                y_tip = offset - 2
                y_text = y_tip - 2
                y_line = y_text - fheight / 2 + fdescent
            else:
                y_tip = offset + 2
                y_text = y_tip + 2 + fheight
                y_line = y_text - fheight / 2 + fdescent
            print y_tip, y_line, y_text
            c.move_to( start_x, y_tip )
            c.line_to( start_x, y_line)
            c.line_to( start_x + (line_length - width)/2 - xbearing -3 , y_line)
            c.move_to( start_x + (line_length - width)/2 - xbearing, y_text)
            c.show_text(label)
            c.move_to( start_x + (line_length + width)/2 + xbearing +3, y_line)
            c.line_to( end_x, y_line)
            c.line_to( end_x, y_tip)
            c.stroke()
        if not lower:
            offset -= (height + 2)
        else:
            offset += (height + 2)

draw_annotation( upper_bars, OFF_Y, c, lower=False )
draw_annotation( lower_bars, OFF_Y + HEIGHT, c)

surface.flush()




