#!/usr/bin/env python3

# Example of use:
# Insert into config template (for plain config file, just replace $ values with actual numbers)
# exec{import generate_brush}
# eval{generate_brush.generate($x0film, $dxbrush, $y0brush, $dybrush, $brushperiod, $y1film, $xnpic_brush, $ynpic_brush, $znpic_brush)}
# Example with horizontal brush:
# eval{generate_brush.generate_horizontal(2, $dxbrush, 2.5, 2.5 + $dybrush, $brushperiod, 13, $xnpic_brush, $ynpic_brush, $znpic_brush)}

def generate(x0film, dxbrush, y0brush, dybrush, brushperiod, ymax, xnpic, ynpic, znpic):
    res = []
    ydispl = 0.0
    while ydispl + y0brush + dybrush <= ymax:
        res.append('film = on')
        res.append('x0film = {0}'.format(x0film - dxbrush))
        res.append('filmwidth = {0}'.format(dxbrush))
        res.append('gradwidth = 0')
        res.append('y0film = {0}'.format(ydispl + y0brush))
        res.append('y1film = {0}'.format(ydispl + y0brush + dybrush))
        res.append('z0film = 0.5')
        res.append('z1film = 10.5')
        res.append('nfilm = 1 ne')
        res.append('mcr = 2.0')
        res.append('xnpic_film = {0}'.format(xnpic))
        res.append('ynpic_film = {0}'.format(ynpic))
        res.append('znpic_film = {0}'.format(xnpic))
        ydispl += brushperiod
    return "\n".join(res)


def generate_horizontal(x0brush, dxbrush, y0brush, y1brush, z0brush, z1brush, brushperiod, xmax, xnpic, ynpic, znpic):
    res = []
    xdispl = 0.0
    while xdispl + x0brush + dxbrush <= xmax:
        res.append('film = on')
        res.append('x0film = {0}'.format(xdispl + x0brush))
        res.append('filmwidth = {0}'.format(dxbrush))
        res.append('gradwidth = 0')
        res.append('y0film = {0}'.format(y0brush))
        res.append('y1film = {0}'.format(y1brush))
        res.append('z0film = {0}'.format(z0brush))
        res.append('z1film = {0}'.format(z1brush))
        res.append('nfilm = 1 ne')
        res.append('mcr = 2.0')
        res.append('xnpic_film = {0}'.format(xnpic))
        res.append('ynpic_film = {0}'.format(ynpic))
        res.append('znpic_film = {0}'.format(xnpic))
        xdispl += brushperiod
    return "\n".join(res)

