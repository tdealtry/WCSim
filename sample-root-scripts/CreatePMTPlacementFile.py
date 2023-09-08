#!/bin/env python3
from pprint import pprint
import math
import ROOT as R

npmts_total = {}
can = R.TCanvas()
drawn = False
leg = R.TLegend()
cols = [R.kBlack, R.kBlue, R.kMagenta]
styles = [31, 20, 27]

def FormatGraph(g, name, pmt_type, title):
    g.SetName(name)
    g.SetTitle()
    g.SetMarkerColor(cols[pmt_type])
    g.SetMarkerStyle(styles[pmt_type])
    g.SetMarkerSize(2)
    leg.AddEntry(g, str(pmt_type), "M")
    
def CalculateCapGrid(npmts, half_height, radius, grid_distance, f, verbose, pmt_type, plot_mode):
    print(f'Attempting to lay out {npmts} PMTs on a cap grid of radius {radius} with ideal PMT distance {grid_distance}')
    #Taking away 50cm here, so as to not overlap with the barrel
    max_radius = radius - 50
    max_radius_sq = max_radius**2
    npmts_new = 0
    x = -max_radius
    while x < +max_radius:
        y = -max_radius
        while y < +max_radius:
            if x**2 + y**2 > max_radius_sq:
                y += grid_distance
                continue
            npmts_new += 1
            y += grid_distance
        x += grid_distance
    frac_diff = (npmts_new - npmts) / npmts
    print(f'{npmts_new} found. This is {frac_diff * 100:.2f}% away from the target number of PMTs')
    #if we're not within x%, try again
    if abs(frac_diff) > 0.035:
        delta = +0.5 if (npmts - npmts_new < 0) else -0.5
        CalculateCapGrid(npmts, half_height, radius, grid_distance + delta, f, verbose, pmt_type, plot_mode)
        return
    #Now we know the grid, lets lay the PMTs out
    g2d = R.TGraph(npmts_new)
    FormatGraph(g2d, f'g2D_CAP_{pmt_type}', pmt_type, ';X;Y')
    g3d = R.TGraph2D(npmts_new * 2)
    FormatGraph(g3d, f'g3D_CAP_{pmt_type}', pmt_type, ';X;Y;Z')
    for icap, (z, cap) in enumerate([[+half_height, 'TOP'], [-half_height, 'BOTTOM']]):
        f.write(f'#DATASTART {cap} CAP NPMTs:{npmts_new} GRID DISTANCE:{grid_distance}\n')
        ipmt = 0
        x = -max_radius
        while x < +max_radius:
            y = -max_radius
            while y < +max_radius:
                if x**2 + y**2 > max_radius_sq:
                    y += grid_distance
                    continue
                if verbose >= 2:
                    print(f'{ipmt} = x:{x:.4f} y:{y:.4f} z:{z:.4f}')
                f.write(f'{x} {y} {z} {0.0} {0.0} {0.0} {pmt_type+1}\n')
                if not icap:
                    g2d.SetPoint(ipmt, x, y)
                g3d.SetPoint(ipmt + icap * npmts_new, x, y, z)
                ipmt += 1
                y += grid_distance
            x += grid_distance
    npmts_total[pmt_type] += npmts_new * 2
    global drawn
    if plot_mode == '2DCap':
        g2d.DrawClone("P" + ("" if drawn else "A"))
        drawn = True
    elif plot_mode == '3D':
        g3d.DrawClone("P" + ("" if drawn else "A"))
        drawn = True
             
def CalculateBarrelGrid(npmts, half_height, radius, grid_distance, f, verbose, pmt_type, plot_mode):
    circumference = 2 * math.pi * radius
    print(f'Attempting to lay out {npmts} PMTs on a barrel grid of circumference {circumference} and total height {half_height * 2} with ideal PMT distance {grid_distance}')
    #Taking away 50cm here, so as to not overlap with the caps
    max_half_height = half_height - 50
    nrows = math.floor((max_half_height - 100) / grid_distance)
    ncols = math.floor(circumference / grid_distance)
    npmts_new = ncols * nrows
    print(f'ncols x nrows = {ncols} x {nrows} = {npmts_new} PMTs using PMT distance {grid_distance}')
    frac_diff = (npmts_new - npmts) / npmts
    print(f'This is {frac_diff * 100:.2f}% away from the target number of PMTs')
    #if we're not within x%, try again
    if abs(frac_diff) > 0.02:
        delta = +0.5 if (npmts - npmts_new < 0) else -0.5
        CalculateBarrelGrid(npmts, half_height, radius, grid_distance + delta, f, verbose, pmt_type, plot_mode)
        return
    #Now we know the grid, lets lay the PMTs out
    f.write(f'#DATASTART BARREL NPMTs:{npmts_new} GRID DISTANCE:{grid_distance}\n')
    g2d = R.TGraph(npmts_new)
    FormatGraph(g2d, f'g2D_BARREL_{pmt_type}', pmt_type, ';PHI;Z')
    g3d = R.TGraph2D(npmts_new)
    FormatGraph(g3d, f'g3D_BARREL_{pmt_type}', pmt_type, ';X;Y;Z')
    ipmt = 0
    for icol in range(ncols):
        phi = icol * (2 * math.pi) / ncols
        x = radius * math.cos(phi)
        y = radius * math.sin(phi)
        for irow in range(nrows):
            #remember the 50 cm offsets, so we don't overlap the caps
            z = - max_half_height + irow * (2 * max_half_height) / nrows
            if verbose >= 2:
                print(f'{icol}, {irow} = phi:{phi:.4f} x:{x:.4f} y:{y:.4f} z:{z:.4f}')
            f.write(f'{x} {y} {z} {0.0} {0.0} {0.0} {pmt_type+1}\n')
            g2d.SetPoint(ipmt, phi, z)
            g3d.SetPoint(ipmt, x, y, z)
            ipmt += 1
    npmts_total[pmt_type] += npmts_new
    global drawn
    if plot_mode == '2DBarrel':
        g2d.DrawClone("P" + ("" if drawn else "A"))
        drawn = True
    elif plot_mode == '3D':
        g3d.DrawClone("P" + ("" if drawn else "A"))
        drawn = True

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='All units in cm')
    parser.add_argument('--npmts-barrel', type=int, nargs='+', default=[12600, 535, 6570/3], help='Target number of PMTs of each type in the barrel region. Recommend to use the order 20", mPMT, OD')
    parser.add_argument('--npmts-cap', type=int, nargs='+', default=[3228, 143, 1508/3], help='Target number of PMTs of each type in a single cap')
    parser.add_argument('--radius', type=float, nargs='+', default=[3248.921111111, 3248.921111111, 3300.927586207], help='Radius to place PMTs on. Taken from a geofile.txt using auto PMT placement')
    parser.add_argument('--half-height', type=float, nargs='+', default=[3296.471111111, 3296.471111111, 3346.477586207], help='Half-Height to place PMTs on. Taken from a geofile.txt using auto PMT placement')
    parser.add_argument('--grid-distance', type=float, nargs='+', default=[100, 350, 200], help='Target grid size for each PMT type')
    parser.add_argument('--filename', type=str, default='PMT_Position.txt', help='Output filename')
    parser.add_argument('--verbose', type=int, default=0, help='Higher number = more text')
    parser.add_argument('--plot-mode', choices=['3D', '2DCap', '2DBarrel'], help='Which plots to plot. Why this option? TCanvas::cd() isn\'t working in my setup...')
    args = parser.parse_args()

    npmt_types = len(args.npmts_cap)
    assert(npmt_types == len(args.npmts_barrel))
    assert(npmt_types == len(args.grid_distance))
    assert(npmt_types == len(args.radius))
    assert(npmt_types == len(args.half_height))
    assert(npmt_types >=1 and npmt_types <= 3)

    print('Attempting to create geometry with:')
    for i in range(3):
        print(f' {2 * args.npmts_cap[i] + args.npmts_barrel[i]} of type {i+1}')

    with open(args.filename, 'w') as f:
        for i in range(npmt_types):
            print('\n\nLooking at PMT type', i)
            npmts_total[i] = 0
            CalculateBarrelGrid(args.npmts_barrel[i], args.half_height[i], args.radius[i], args.grid_distance[i], f, args.verbose, i, args.plot_mode)
            CalculateCapGrid(args.npmts_cap[i], args.half_height[i], args.radius[i], args.grid_distance[i], f, args.verbose, i, args.plot_mode)

    print('\nTotal number of PMTs:')
    pprint(npmts_total)
    if args.plot_mode:
        can.BuildLegend()
