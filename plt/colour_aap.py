#!/usr/bin/env python

import sys
import math
from bioio import read_fasta, MultiSeqs

# TODO:
# set antialias when writing text! (meaning we need to draw everything first,
# and then write text after that)

LINEWIDTH = 11
LINEHEIGHT = 14
OFFSET = 120
OFF_Y = 100
PAD_BTM = 50

aa_colors = {
    'H': ( 155, 221, 255),  # columbia blue
    'K': ( 0, 112, 255),  # brandeis blue
    'R': ( 0, 255, 255),  # aqua
    'D': ( 255, 0, 0),      # red
    'E': ( 220, 20, 60),    # crimson
    'F': ( 238, 130, 238),  # violet (alt: indigo)
    'W': ( 0, 0, 255),      # blue
    'Y': ( 65, 105, 225),   # royal blue (webcolor)
    'T': ( 156, 186, 227),  # carolina blue
    'S': ( 251, 174, 210),  # lavender pink
    'N': ( 223, 115, 255),  # hellotrope
    'Q': ( 128, 0, 128),    # purple
    'I': ( 34, 139, 34),    # forest green
    'L': ( 107, 142, 35),        # olive drab
    'V': ( 154, 205, 50),   # yellow-green
    'A': ( 0, 255, 127),    # spring green
    'G': ( 0, 255, 0),      # green (or harlequin)
    'P': ( 255, 127, 0),    # orange
    'M': ( 112, 66, 20),    # sepia
    'C': ( 255, 255, 0),    # yellow
    'X': ( 255, 255, 255)   # white
}

aa_aesnn3 = {
    'A': ( -0.99, -0.61,  0.00 ),
    'R': (  0.28, -0.99, -0.22 ),
    'N': (  0.77, -0.24,  0.59 ),
    'D': (  0.74, -0.72, -0.35 ),
    'C': (  0.34,  0.88,  0.35 ),
    'Q': (  0.12, -0.99, -0.99 ),
    'E': (  0.59, -0.55, -0.99 ),
    'G': ( -0.79, -0.99,  0.10 ),
    'H': (  0.08, -0.71,  0.68 ),
    'I': ( -0.77,  0.67, -0.37 ),
    'L': ( -0.92,  0.31, -0.99 ),
    'K': ( -0.63,  0.25,  0.50 ),
    'M': ( -0.80,  0.44, -0.71 ),
    'F': (  0.87,  0.65, -0.53 ),
    'P': ( -0.99, -0.99, -0.99 ),
    'S': (  0.99,  0.40,  0.37 ),
    'T': (  0.42,  0.21,  0.97 ),
    'W': ( -0.13,  0.77, -0.90 ),
    'Y': (  0.59,  0.33, -0.99 ),
    'V': ( -0.99,  0.27, -0.52 ),
    'X': (  0.99,  0.99,  0.99 )
}

white = (1, 1, 1)
black = (0, 0, 0)

aa_color_fg = {
    'A': white,
    'R': white,
    'N': black,
    'D': black,
    'C': black,
    'Q': white,
    'E': white,
    'G': white,
    'H': black,
    'I': black,
    'L': white,
    'K': black,
    'M': white,
    'F': black,
    'P': white,
    'S': black,
    'T': black,
    'W': black,
    'Y': black,
    'V': white,
    'X': black
}
    

def normalize_aesnn3_color():
    global aa_colors, aa_aesnn3
    for k in aa_aesnn3:
        r, g, b = aa_aesnn3[k]
        aa_colors[k] = ( (r+1)/2, (g+1)/2, (b+1)/2 )


def normalize_color():
    global aa_colors
    for k in aa_colors:
        r, g, b = aa_colors[k]
        aa_colors[k] = ( 1.0 * r / 255, 1.0 * g / 255, 1.0 * b / 255 )

def create_spectrum(reg_size):

    step = math.pi/(reg_size - 1)
    #colours = [ ( (-math.cos(x*step)), (1-math.cos(x*step*6))/2, (math.cos(x*step)) ) for x in range(reg_size) ]
    #colours = [ (r*255, g*255, b*255) for (r,g,b) in colours ]
    colours = [ ( -math.cos(x*step)*1.25, (1-math.cos(x*step*6))/2, math.cos(x*step)*1.25 ) for x in range(reg_size) ]
    return colours


def draw_label( offset, labels ):
    pass


def draw_aamap(mseq, m_len, s_len, outfile, upper_label=None, reg_pos=None, tolerance=None):

    import cairo

    regions = []
    for i in reg_pos:
        (regname, regpos) = i.split(':')
        regpos = int(regpos)
        regions.append( (regname, regpos) )

    W, H = s_len * LINEWIDTH + OFFSET, m_len * LINEHEIGHT + OFF_Y + PAD_BTM
    #surface = cairo.ImageSurface( cairo.FORMAT_ARGB32, W, H )
    surface = cairo.PSSurface( outfile, W, H )

    c = cairo.Context(surface)

    # clear background
    c.set_antialias( cairo.ANTIALIAS_NONE )
    c.set_line_cap(cairo.LINE_CAP_BUTT)
    c.set_source_rgb( 1, 1, 1 )
    c.set_operator(cairo.OPERATOR_SOURCE)
    c.paint()

    # draw upper label

    #off_y = OFF_Y - 5

    #c.set_source_rgb(0, 0, 0)
    #c.set_font_size(LINEWIDTH-1)
    #fascent, fdescent, fheight, fxadvance, fyadvance = c.font_extents()
    #for pos, label in enumerate(upper_label):
    #    print pos, label
    #    c.move_to( OFFSET + LINEWIDTH*(pos) + 3, off_y )
    #    c.save()
    #    c.rotate( - math.pi / 2)
    #    c.show_text('%d' % label)
    #    c.restore()

    # draw aa region

    off_y = OFF_Y

    for y in range(m_len):

        c.set_line_width(LINEWIDTH)
        s = mseq[y]

        c.set_font_size(LINEWIDTH-3)
        for x in range(s_len):

            c.move_to(OFFSET + LINEWIDTH*x, off_y)
            aa = s.get_seq()[x]
            c.set_source_rgb( *aa_colors[aa] )
            c.line_to(OFFSET + LINEWIDTH*x, off_y + LINEHEIGHT)
            c.stroke()
            c.set_source_rgb( *aa_color_fg[aa] )
            fascent, fdescent, fheight, fxadvance, fyadvance = c.font_extents()
            xbearing, ybearing, w, h, xadv, yadv = c.text_extents( aa )
            c.move_to(OFFSET + LINEWIDTH * x - w/2 + xbearing - 1, off_y + (LINEHEIGHT + fheight)/2 - 1 )
            c.show_text( aa )
            c.stroke()
            

        c.set_font_size(LINEHEIGHT-2)
        c.set_source_rgb(0, 0, 0)
        fascent, fdescent, fheight, fxadvance, fyadvance = c.font_extents()
        c.move_to( 5, off_y + (LINEHEIGHT + fheight)/2 - 1 )
        c.show_text( s.name.replace('_',' ') )

        off_y += LINEHEIGHT

    #surface.write_to_eps(outfile)
    #surface.paint()

    # draw tolerance
    if tolerance:
        off_y += 10

        c.move_to( 5, off_y + (LINEHEIGHT + fheight)/2 - 1)
        c.show_text( 'Tolerance' )

        #c.set_line_width(1)
        #c.move_to(OFFSET, off_y + LINEHEIGHT/2)
        #c.line_to(OFFSET + LINEWIDTH*s_len, off_y + LINEHEIGHT/2)
        #c.stroke()

        c.set_line_width(LINEWIDTH)
        for x in range(s_len):
            c.move_to(OFFSET + LINEWIDTH*x, off_y + LINEHEIGHT/2)
            if tolerance[x] == 0:
                next_y = off_y + LINEHEIGHT/2 + 5
                c.set_source_rgb(255, 0, 0)
            else:
                next_y = off_y + LINEHEIGHT/2 - 5
                c.set_source_rgb(0, 255, 0)
            c.line_to(OFFSET + LINEWIDTH*x, next_y)
            c.stroke()

    #off_y += 25
    off_y = OFF_Y - 35

    # draw region of protein

    new_points = []
    M = len(regions)
    start_x = OFFSET - LINEWIDTH/2
    c.set_line_width( LINEHEIGHT )
    r,g,b = 0, 0, 1
    #c.set_font_size(LINEWIDTH - 1)
    c.set_font_size(LINEWIDTH + 6)
    fascent, fdescent, fheight, fxadvance, fyadvance = c.font_extents()
    reg_colours = create_spectrum(M)
    for x in range(M):
        c.set_source_rgb(*reg_colours[x])
        #c.set_line_width( LINEHEIGHT )
        regname, startpos = regions[x]
        if x+1 < M:
            endpos = regions[x+1][1]
        else:
            endpos = upper_label[-1]
        end_point = 0
        for p in upper_label:
            if p < startpos: continue
            if p > endpos: break
            end_point += 1
            new_points.append( p - startpos + 1 )
        if end_point == 0: continue
        c.move_to( start_x, off_y )
        end_x = start_x + end_point * LINEWIDTH
        c.line_to( end_x, off_y )
        text_point = (end_x - start_x)/2 + start_x + 0.5 * fheight
        c.move_to( text_point, off_y - 12 )
        c.save()
        c.rotate( - math.pi / 2 )
        c.set_source_rgb(0,0,0)
        c.show_text(regname)
        c.restore()
        c.stroke()
        r,g,b = r + 0.1, g + 0.025, b - 0.1
        start_x = end_x

    off_y = OFF_Y - 5

    c.set_source_rgb(0, 0, 0)
    c.set_font_size(LINEWIDTH + 1)
    fascent, fdescent, fheight, fxadvance, fyadvance = c.font_extents()
    for pos, label in enumerate(new_points):
        #print pos, label
        c.move_to( OFFSET + LINEWIDTH*(pos) + 3, off_y )
        c.save()
        c.rotate( - math.pi / 2)
        c.show_text('%d' % label)
        c.restore()

    surface.flush()
    #surface.write_to_png(sys.argv[2])


def usage():
    print "Usage: colour_aap.py infile_fas outfile_eps"
    print
    sys.exit(0)

if __name__ == '__main__':

    if len(sys.argv) < 2:
        usage()

    outfile = sys.argv[2]
    mseq = read_fasta(sys.argv[1], MultiSeqs( reserved_names = [ '_tag' ]))

    s_len = 0
    inseqs = MultiSeqs()
    upper_label = None
    tolerance = None
    regpos = []
    for s in mseq:
        if s.name.startswith('_tag'):
            if 'position' in s.get_attribute():
                upper_label = [ int(x) for x in s.seq.split(',') ]
            elif 'region' in s.get_attribute():
                regpos = s.seq.split(',')
            elif 'tolerance' in s.get_attribute():
                tolerance = [ int(x) for x in s.seq ]
            continue
        inseqs.append( s )

    normalize_aesnn3_color()
    draw_aamap( inseqs, len(inseqs), len(inseqs[0]), outfile, upper_label, regpos, tolerance)

        
