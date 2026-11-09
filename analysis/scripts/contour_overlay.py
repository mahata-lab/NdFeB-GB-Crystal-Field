"""Fixed physical-density contour overlay for the current Figure 5."""
import numpy as np
from PIL import Image
from reportlab.lib.colors import HexColor

CONTOUR_LEVELS=(0.25,0.75,2.0)  # same Å^-2 values in all six panels

def contours(c, density, x, y, width, height):
    """Marching-squares isolines at fixed physical-density values on the displayed grid."""
    ny, nx = density.shape
    c.setStrokeColor(HexColor('#202020'))
    c.setLineWidth(0.32)
    for level in CONTOUR_LEVELS:
        bits = ((density[:-1,:-1] >= level).astype(np.uint8)
                + 2*(density[:-1,1:] >= level)
                + 4*(density[1:,1:] >= level)
                + 8*(density[1:,:-1] >= level))
        yy, xx = np.nonzero((bits != 0) & (bits != 15))
        path = c.beginPath()
        for j, i in zip(yy, xx):
            values = [density[j,i], density[j,i+1], density[j+1,i+1], density[j+1,i]]
            vertices = [(i,j),(i+1,j),(i+1,j+1),(i,j+1)]
            points = []
            for k in range(4):
                n = (k+1) % 4
                if (values[k] >= level) == (values[n] >= level):
                    continue
                ratio = (level-values[k])/(values[n]-values[k])
                px = vertices[k][0] + ratio*(vertices[n][0]-vertices[k][0])
                py = vertices[k][1] + ratio*(vertices[n][1]-vertices[k][1])
                points.append((x+px*width/(nx-1), y+py*height/(ny-1)))
            if len(points) == 2:
                pairs = [(0,1)]
            elif len(points) == 4:
                # Resolve a saddle by the bilinear cell-centre value.
                pairs = [(0,1),(2,3)] if (values[0]>=level)==(sum(values)/4>=level) else [(0,3),(1,2)]
            else:
                continue
            for a,b in pairs:
                path.moveTo(*points[a]); path.lineTo(*points[b])
        c.drawPath(path)


def draw_contours(c,density,x,y,width,height):
    # Refine only the displayed field; all numerical integrals use the original 3D cubes.
    smooth=np.maximum(np.asarray(Image.fromarray(density.astype('float32')).resize((220,220),Image.Resampling.BICUBIC)),0)
    contours(c,smooth,x,y,width,height)
