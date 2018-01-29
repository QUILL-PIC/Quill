#!/usr/bin/python

from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import resread
import tcmap
from collections import OrderedDict
import expression_parser


def __get_data_folder(data_folder, kwargs_dict):
    if data_folder is not None:
        return data_folder
    elif 'df' in kwargs_dict:
        return kwargs_dict.pop('df')
    else:
        return None


def tex_format(space_item):
    tmp = space_item
    if space_item in ['ex','ey','ez','bx','by','bz']:
        tmp = space_item[0].upper() + '_' + space_item[1]
    elif space_item in ['ux','uy','uz']:
        tmp = 'p_' + space_item[1]
    elif space_item in ['vx','vy','vz']:
        tmp = space_item[0] + '_' + space_item[1]
    elif space_item in ['xi','chi','theta']:
        tmp = '\\' + space_item
    elif space_item == 'g':
        tmp = '\\gamma'
    elif space_item == 'phi':
        tmp = '\\varphi'
    return '$' + tmp + '$'


def density(t=0, plane='xy', max_w=0, max_e_density=0, max_p_density=0, max_g_density=0, max_i_density=0, axis=None,
            extent=None, save2=None, data_folder=None, particles='geipw', cmaps={}, xlim=None, ylim=None, clf=False,
            color_mixing = False, **kwargs):
    """
    Plots density distributions of particles and the electromagnetic field w.

    Parameters
    ----------
    plane : {'xy', 'xz', 'yz'}
    axis
        is passed to plt.axis method.
    extent
        Passed to imshow if not None or 'xi'. If None, extent is calculated automatically.
        If 'xi', will use xi = x - vt instead of x.
    save2 : string
        File to save the image to. If None, does nothing.
    data_folder : string
        E.g. '../results/'. Use 'df' for shortcut.
    particles : string
        String of particles to plot (e.g. 'gie', 'we'). Possible values are: g, e, i, p, w.
        The order determines the plot order (e.g. 'ew' and 'we' have different plot orders).
    cmaps : dict
        Dictionary of particle-colormap pairs to be used instead of defaults. E.g. {'e': 'Blues', 'w': 'tcmap_green'}.
    color_mixing : bool
        Turns on color mixing instead of overlay of palettes for multiple particles. The default palettes are changed to
        opaque ones.
    kwargs
        parameters to be passed to the underlying 'imshow' methods.
    """
    if color_mixing:
        import qcolor

    resread.t = '%g' % t
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()

    if clf:
        plt.clf()
    plt.gca().set_title('Particle densities', fontsize='medium')
    plt.xlabel(tex_format(plane[0])[:-1] + '/\lambda$')
    plt.ylabel(tex_format(plane[1])[:-1] + '/\lambda$')
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if plane == 'xz':
        xlength = resread.nx * resread.dx
        ylength = resread.nz * resread.dz
    elif plane == 'yz':
        xlength = resread.ny * resread.dy
        ylength = resread.nz * resread.dz
    else:
        xlength = resread.nx * resread.dx
        ylength = resread.ny * resread.dy
    if axis is not None:
        plt.axis(axis)

    if 'e' in particles:
        edensity = - resread.density('rho', plane)
        if max_e_density == 0:
            max_e_density = edensity.max()
            print('qplot.density: max_e_density = {0}'.format(max_e_density))

    if 'w' in particles:
        w = resread.density('w', plane)
        if max_w == 0:
            max_w = w.max()
            print('qplot.density: max_w = {0}'.format(max_w))

    if max_p_density == 0:
        max_p_density = max_e_density
    if max_g_density == 0:
        max_g_density = max_e_density
    if max_i_density == 0:
        max_i_density = max_e_density

    if extent is None:
        extent = (0, xlength, 0, ylength)
    elif extent == 'xi':
        extent = (resread.xi(0, t), resread.xi(xlength, t), 0, ylength)

    tmp_kwargs = {'extent': extent, 'vmin': 0,  'origin': 'lower', 'interpolation': 'none'}
    tmp_kwargs.update(kwargs)
    kwargs = tmp_kwargs

    if color_mixing:
        tmp_cmaps = {'g': 'qcolor_blue', 'w': 'qcolor_orange', 'i': 'qcolor_purple', 'e': 'qcolor_green', 'p': 'qcolor_red'}
    else:
        tmp_cmaps = {'g': 'tcmap_blue', 'w': 'tcmap_orange', 'i': 'tcmap_purple', 'e': 'tcmap_green', 'p': 'tcmap_red'}
    tmp_cmaps.update(cmaps)
    cmaps = tmp_cmaps

    for p in particles:
        if p == 'g' and 'g' in resread.particles_for_output:
            g_density = resread.density('rho_ph', plane)
            if max_g_density == 0:
                max_g_density = g_density.max()
                print('qplot.density: max_g_density = {0}'.format(max_g_density))
            plt.imshow(g_density, cmap=cmaps['g'], vmax=max_g_density, **kwargs)
        elif p == 'w':
            plt.imshow(w, cmap=cmaps['w'], vmax=max_w, **kwargs)
        elif p == 'i' and resread.icmr != []:
            i_density = resread.density('irho_' + str(resread.icmr[0]) + '_', plane)
            if max_i_density == 0:
                max_i_density = i_density.max();
                print('qplot.density: max_i_density = {0}'.format(max_i_density))
            plt.imshow(i_density, cmap=cmaps['i'],
                       vmax=max_i_density, **kwargs)
        elif p == 'e' and 'e' in resread.particles_for_output:
            plt.imshow(edensity, cmap=cmaps['e'], vmax=max_e_density, **kwargs)
        elif p == 'p' and 'p' in resread.particles_for_output:
            p_density = resread.density('rho_p', plane)
            if max_p_density == 0:
                max_p_density = p_density.max()
                print('qplot.density: max_p_density = {0}'.format(max_p_density))
            plt.imshow(p_density, cmap=cmaps['p'], vmax=max_p_density, **kwargs)
        else:
            if p not in 'gwiep':
                print("Warning: Particle '%c' is not a valid option. Possible options are: gwiep." % p)

    if color_mixing:
        imgs = [x for x in plt.gca().get_children() if isinstance(x, mpl.image.AxesImage)]
        rgb_values = qcolor.mix_images(imgs)
        for img in imgs:
            img.remove()
        plt.imshow(rgb_values, **kwargs)

    if save2 is not None:
        plt.savefig(save2)


def filter_phspace(phs, expr, var_space):
    mask = np.ones(phs.shape[1], dtype=bool).T
    if expr:
        mask &= expression_parser.evaluate(expr, lambda var_name: phs.T[:, var_space.index(var_name)])
    return phs.T[mask].T


def space_to_list(space, filter_expr):
    params_list = ['ux', 'uy', 'uz', 'vx', 'vy', 'vz', 'theta', 'phi', 'chi', 'xi', 'x', 'y', 'z', 't', 'g', 'q']
    if isinstance(space, str):
        space_tmp = space
        space_dict = OrderedDict({})
        for token in params_list:
            if space_tmp.find(token) > -1:
                space_dict[space_tmp.find(token)] = token
                space_tmp = space_tmp.replace(token, ' ' * len(token)) # need to preserve params absolute positions in space_tmp
        if len(space_tmp.replace(' ', '')) > 0:
            raise ValueError('Invalid space specified: {0}'.format(space))
        space = list(OrderedDict(sorted(space_dict.items())).values())
    return list(OrderedDict((v, i) for i, v in enumerate(space)).keys())


def particles(t=0, space=['x','y'], particles='geip', colors='bgmrcyk', r=3, alpha=0.1, cmap='jet', gamma=0,
              data_folder=None, axis=[], save2=None, vmin=None, vmax=None, xlim=None, ylim=None, filter=None, clf=False,
              include_deleted=False, every=1, **kwargs):
    """
    Plots particles as dots in (phase)*space*.
    Examples:
        qplot.particles(10,[\'x\',\'y\']),\n\
        qplot.particles(15,\'xyg\',\'i\',gamma=1).
    :param t: time for which to plot the particle distribution
    :param space: phasespace containing 2 or 3 variables (if 3, the third one is shown as color on a 2D plot). Can be either list or string
    :param particles: particles to plot (g for photons, e for electrons, i for ions
    :param filter: condition for filtering the variables. E.g.: 'x > 5 && (theta > 0.5 || theta < -0.5)'
        Allowed operations: + - * / ^ () < > = != <= >= && || sin cos exp log min max
            (there are also aliases '**','and','or')
    """

    def particle_from_suffix(s):
        if s == '':
            return 'e'
        elif s == '_p':
            return 'p'
        elif s == '_ph':
            return 'g'
        else:
            return 'i'

    resread.t = '%g' % t
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()

    filter_expr = expression_parser.to_polish(filter)
    space = space_to_list(space, filter_expr)
    if len(space) == 2:
        is_3d = False
    elif len(space) == 3:
        is_3d = True
    else:
        raise ValueError('Incorrect value for space parameter')
    space += expression_parser.get_vars(filter_expr)

    if clf:
        plt.clf()
    plt.xlabel(tex_format(space[0]))
    plt.ylabel(tex_format(space[1]))
    if axis:
        plt.axis(axis)
    lp = list(particles)
    s = []
    for p in lp:
        if p == 'e':
            s.append('')
        elif p == 'p':
            s.append('_p')
        elif p == 'g':
            s.append('_ph')
        elif p == 'i':
            for cmr in resread.icmr:
                s.append('_'+str(cmr)+'_')
    if not is_3d:
        c = list(colors)
        i = 0
        for suffix in s:
            phspace = resread.particles('phasespace'+suffix, space, every=every)
            if (resread.catching or (resread.dump_photons and suffix == '_ph')) and include_deleted:
                for t1 in np.arange(0, t+0.001, resread.output_period):
                    resread.t = '%g' % t1
                    phspace_del = resread.particles('deleted'+suffix, space, every=every)
                    phspace = np.append(phspace, phspace_del, axis=1)
                resread.t = '%g' % t
            phspace = filter_phspace(phspace, filter_expr, space)
            print('Plotting {0}({1}) for {2}; vmin = {3}, vmax = {4}'
                  .format(space[1], space[0], particle_from_suffix(suffix), vmin, vmax))
            plt.scatter(phspace[0,:], phspace[1,:], color=c[i], s=r*r, alpha=alpha, edgecolors='None')
            if ylim:  # we assume that either ylim is set or vmin / vmax are set
                plt.ylim(ylim)
            else:
                plt.ylim(vmin,vmax)
            plt.xlim(xlim)
            i += 1
    else:
        cmap = tcmap.get(cmap, gamma=gamma)
        plt.title(tex_format(space[2]))
        phspace = resread.particles('phasespace'+s[0], space, every=every)  # 2 or more sorts of ions are not supported
        if (resread.catching or (resread.dump_photons and s[0] == '_ph')) and include_deleted:
            for t1 in np.arange(0, t+0.001, resread.output_period):
                resread.t = '%g' % t1
                phspace_del = resread.particles('deleted'+s[0], space, every=every)
                phspace = np.append(phspace, phspace_del, axis=1)
            resread.t = '%g' % t
        phspace = filter_phspace(phspace, filter_expr, space)
        if vmin is None:
            vmin = np.min(phspace[2,:])
        if vmax is None:
            vmax = np.max(phspace[2,:])
        print('Plotting {0}({1},{2}) for {3}; vmin = {4}, vmax = {5}'
              .format(space[2], space[0], space[1], particles[0], vmin, vmax))
        plt.scatter(phspace[0,:], phspace[1,:], s=r*r, c=phspace[2,:],
                    cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='None')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.colorbar()
    if save2 is not None:
        plt.savefig(save2)


def dist_fn(t=0, data_folder=None, particles='e', space='x', energy=False, nbins=200, filter=None, clf=False, global_limits=True,
            log_scale=None, save2=None, every=1, include_deleted=False, **kwargs):
    """
    Plots N and energy distribution functions over the specified space.

    Parameters
    ----------
    data_folder : string
        E.g. '../results/'. Use 'df' for shortcut.
    particles : string
        String of particle species to plot: 'e', 'p', 'g', or 'i'. Only one particle species at time is allowed.
    space : string, list
        One- or two-dimensional space over which the distribution function is calculated.
    energy: bool
        if False (by default), N (number of particles) distibution over space is shown.
        if True, energy distribution is shown instead
    nbins: int
        The number of bins for calculating the histogram.
    filter: string
        Expression to filter the particles by. E. g.: 'x > 5 && y < 10', 'g > 0.1*max(g)'
    global_limits: bool
        If true, the distribution range is set to global min / max for each parameter (not for min / max of particles)
        E. g. x ranges from 0 to nx*dx, phi --- from -pi to pi.
    log_scale: string
        None, 'x', 'y' or 'xy'
    save2 : string
        File to save the image to. If None, does nothing.
    every: int
        Tells resread to take every n-th file (1 by default)
    kwargs
        parameters to be passed to the plotting methods (hist, imshow).
    """
    resread.t = '%g' % t
    if len(particles) != 1:
        raise ValueError('dist_fn supports only one type of particles')
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()

    if particles == 'i' and len(resread.icmr) == 0:
        raise ValueError('Ions data are missing')
    file_suffixes = {'e': '', 'p': '_p', 'g': '_ph', 'i': '_' + str(resread.icmr[0]) + '_'}
    caption_suffixes = {'e' : '_e', 'p': '_p', 'g': '_{ph}', 'i': '_i'}
    c_suffix = caption_suffixes[particles[0]]
    filter_expr = expression_parser.to_polish(filter)
    space = space_to_list(space, filter_expr)
    is_2d = len(space) > 1
    space.append('q')  # need also to fetch pseudoparticle weight
    if energy:
        space.append('g')
    space += expression_parser.get_vars(filter_expr)
    phspace = resread.particles('phasespace' + file_suffixes[particles[0]], space, every=every)
    if (resread.catching or (resread.dump_photons and particles[0] == 'g')) and include_deleted:
        for t1 in np.arange(0, t+0.001, resread.output_period):
            resread.t = '%g' % t1
            phspace_del = resread.particles('deleted' + file_suffixes[particles[0]], space, every=every)
            phspace = np.append(phspace, phspace_del, axis=1)
        resread.t = '%g' % t
    phspace = filter_phspace(phspace, filter_expr, space)

    if clf:
        plt.clf()
    limits_per_space = {
        'x' : [0, resread.nx * resread.dx], 'y' : [0, resread.ny * resread.dy], 'z' : [0, resread.nz * resread.dz],
        'phi' : [-np.pi, np.pi], 'theta' : [0, np.pi]
    }
    if global_limits:
        if is_2d:
            hist_range = [limits_per_space[space[0]] if space[0] in limits_per_space else [np.min(phspace[0,:]), np.max(phspace[0,:])], 
                limits_per_space[space[1]] if space[1] in limits_per_space else [np.min(phspace[1,:]), np.max(phspace[1,:])]]
        else:
            hist_range = limits_per_space[space[0]] if space[0] in limits_per_space else None
    else:
        hist_range = None

    if is_2d:
        weight = np.abs(phspace[2,:]) * (phspace[3,:] - (1.0 if particles in 'ep' else 0.0) if energy else 1.0)
        hist, xedges, yedges = np.histogram2d(phspace[0,:], phspace[1,:], range=hist_range, bins=nbins, weights=weight, normed=True)
        plt.imshow(hist.T, origin='lower', aspect='auto', extent=(xedges[1], xedges[-1], yedges[1], yedges[-1]), **kwargs)
        plt.colorbar()
        plt.ylabel(tex_format(space[1]))
        plt.title(tex_format(r'\varepsilon'+c_suffix if energy else 'N'+c_suffix) + '(' + tex_format(space[0]) + ', ' + tex_format(space[1]) + ')')
    else:
        weight = np.abs(phspace[1,:]) * (phspace[2,:] if energy else 1.0)
        #hist, bin_edges = np.histogram(phspace[0,:], range=hist_range, bins=nbins, weights=weight, density=True)
        #plt.plot(bin_edges[1:], hist)
        plt.hist(phspace[0,:], range=hist_range, bins=nbins, weights=weight, normed=True, **kwargs)
        plt.ylabel(tex_format(r'\varepsilon'+c_suffix if energy else 'N'+c_suffix) + '(' + tex_format(space[0]) + ')')
    plt.xlabel(tex_format(space[0]))
    if not log_scale is None and 'x' in log_scale:
        plt.xscale('log')
    if not log_scale is None and 'y' in log_scale:
        plt.yscale('log')
    if save2 is not None:
        plt.savefig(save2)


def tracks(space=['x','y'], particles='geip', t0=0, t1=0, colors='bgmrcyk', cmaps=['jet'], clims2all=1, axis=[], save2=None,
            r=2, data_folder=None, clf=False, every=1, **kwargs):
    'Plots particle tracks as lines in 2D or dots in 3D (phase)*space*\n\
    at [tr_start+*t0*,tr_start+*t1*]\n\
    Examples:\n\
    tracks(space=[\'t\',\'x\',\'g\'],particles=\'e\',clims2all=1)\n\
    tracks(colors=\'\',axis=[9,11,8,12],save2=\'/tmp/tmp.png\') # plots\
    each track with its own color and saves picture to a file\n\
    '
    # clims2all=1 -> color limits for cmaps[i] (vmin,vmax) the same
    # for all particles (trajectories) of this type
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    track_names = []
    for filename in os.listdir(resread.data_folder):
        if filename.find('track')==0:
            track_names.append(filename)
    tracks = []
    cmr = []
    for trackname in track_names:
        resread.t = ''
        tmp = resread.particles(trackname, space, every=every)
        if t1!=0 and int(np.floor(t1/resread.dt)) < len(tmp[0,:]):
            tracks.append(tmp[:,int(np.floor(t0/resread.dt)):int(np.floor(t1/resread.dt))])
        else:
            tracks.append(tmp[:,int(np.floor(t0/resread.dt)):])
        cmr.append(float( trackname[ trackname.find('_')+1 : trackname.find('_',trackname.find('_')+1) ] ))
    def cmr2int(a):
        'Converts cmr *a* to the int index of *colors* or *cmaps*'
        if len(space)==2:
            if colors=='':
                return 0
            else:
                maxv = len(colors)
        elif len(space)==3:
            maxv = len(cmaps)
        if particles.find('i')!=-1:
            s = particles[:particles.find('i')] + (len(resread.icmr)-1)*'i' + particles[particles.find('i'):]
        else:
            s = particles
        if a==-1 and s.find('e')!=-1:
            return s.find('e')%maxv
        elif a==1 and s.find('p')!=-1:
            return s.find('p')%maxv
        elif a==0 and s.find('g')!=-1:
            return s.find('g')%maxv
        elif s.find('i')!=-1:
            for j in np.arange(len(resread.icmr)):
                if a==resread.icmr[j]:
                    return (particles.find('i')+j)%maxv
            return 0
        else:
            return 0

    if clf:
        plt.clf()
    if len(space)==2:
        if colors!='':
            for i in np.arange(len(track_names)):
                if cmr[i]==-1 and particles.find('e')!=-1:
                    plt.plot(tracks[i][0,:],tracks[i][1,:],color=colors[cmr2int(cmr[i])])
                elif cmr[i]==1 and particles.find('p')!=-1:
                    plt.plot(tracks[i][0,:],tracks[i][1,:],color=colors[cmr2int(cmr[i])])
                elif cmr[i]==0 and particles.find('g')!=-1:
                    plt.plot(tracks[i][0,:],tracks[i][1,:],color=colors[cmr2int(cmr[i])])
                elif particles.find('i')!=-1:
                    for j in np.arange(len(resread.icmr)):
                        if cmr[i]==resread.icmr[j]:
                            plt.plot(tracks[i][0,:],tracks[i][1,:],color=colors[cmr2int(cmr[i])])
        else:
            for i in np.arange(len(track_names)):
                if cmr[i]==-1 and particles.find('e')!=-1:
                    plt.plot(tracks[i][0,:],tracks[i][1,:])
                elif cmr[i]==1 and particles.find('p')!=-1:
                    plt.plot(tracks[i][0,:],tracks[i][1,:])
                elif cmr[i]==0 and particles.find('g')!=-1:
                    plt.plot(tracks[i][0,:],tracks[i][1,:])
                elif particles.find('i')!=-1:
                    for j in np.arange(len(resread.icmr)):
                        if cmr[i]==resread.icmr[j]:
                            plt.plot(tracks[i][0,:],tracks[i][1,:])
    elif len(space)==3:
        plt.title(tex_format(space[2]))
        if clims2all==1:
            def pname(a):
                if a==-1:
                    return 'e'
                elif a==1:
                    return 'p'
                elif a==0:
                    return 'g'
                else:
                    return 'i'
            for m in [-1,1,0]+resread.icmr:
                if particles.find(pname(m))!=-1:
                    tmp0 = []
                    tmp1 = []
                    tmp2 = []
                    for i in np.arange(len(track_names)):
                        if cmr[i]==m:
                            for j in np.arange(len(tracks[i][0,:])):
                                tmp0.append(tracks[i][0,j])
                                tmp1.append(tracks[i][1,j])
                                tmp2.append(tracks[i][2,j])
                    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(m)],s=r*r,edgecolors='none')
        else:
            for i in np.arange(len(track_names)):
                tmp0 = []
                tmp1 = []
                tmp2 = []
                for j in np.arange(len(tracks[i][0,:])):
                    tmp0.append(tracks[i][0,j])
                    tmp1.append(tracks[i][1,j])
                    tmp2.append(tracks[i][2,j])
                if particles.find('e')!=-1 and cmr[i]==-1:
                    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,edgecolors='none')
                if particles.find('p')!=-1 and cmr[i]==1:
                    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,edgecolors='none')
                if particles.find('g')!=-1 and cmr[i]==0:
                    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,edgecolors='none')
                if particles.find('i')!=-1:
                    for k in np.arange(len(resread.icmr)):
                        if cmr[i]==resread.icmr[k]:
                            plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,edgecolors='none')
    else:
        print('qplot.tracks: warning: ambiguous value for *space*')
    if axis!=[]:
        plt.axis(axis)
    plt.xlabel(tex_format(space[0]))
    plt.ylabel(tex_format(space[1]))
    if save2 is not None:
        plt.savefig(save2)


def rpattern(t=None, particles='geip', colors='bgmrcyk', dphi=0.1, save2=None, data_folder=None, catching=True, polar=True, clf=False, **kwargs):
    'Plots radiation pattern of the emitted energy\n\
    Examples:\n\
    rpattern() # plots radiation patterns for all particles\n\
    at t_end'
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    if t==None:
        t = resread.t_end - resread.dt
    resread.t = '%g' % t
    lp = list(particles)
    s = []
    suffix_mapping = {'':'e', '_p':'p', '_ph':'g'}
    for p in lp:
        if p=='e':
            s.append('')
        elif p=='p':
            s.append('_p')
        elif p=='g':
            s.append('_ph')
        elif p=='i':
            for cmr in resread.icmr:
                s.append('_'+str(cmr)+'_')
                suffix_mapping[s[-1]] = 'i'
    c = list(colors)
    ci = 0
    if polar:
        plt.subplot(111, polar=True)
    for suffix in s:
        print ('Building rpattern for ' + suffix_mapping[suffix] + ' ...')
        qgphi = resread.particles('phasespace'+suffix,['q','g','phi'])
        phi = np.arange(-np.pi,np.pi+dphi,dphi)
        n = len(phi)
        rp = np.zeros(n)
        for i in np.arange(len(qgphi[0,:])):
            j = int(np.floor((qgphi[2,i]+np.pi)/dphi))
            rp[j] += np.fabs(qgphi[0,i])*qgphi[1,i]
        
        # Including particles that have been deleted at boundaries
        rp_withcatching = None
        if (resread.catching or (resread.dump_photons and suffix == '_ph')) and include_deleted:
            rp_withcatching = np.copy(rp)
            print ('Processing files with deleted particles')
            t_files = ['%g'%t1 for t1 in np.arange(0, t+0.0001, resread.output_period)]
            for t_file in t_files:
                resread.t = t_file
                qgphi_del = np.zeros((3,1))
                try:
                    qgphi_del = resread.particles('deleted' + suffix, ['q','g','phi'])
                except:
                    print ('Unable to access the file [deleted{0}{1}]'.format(suffix, t_file))
                for i in np.arange(len(qgphi_del[0,:])):
                    j = int(np.floor((qgphi_del[2,i]+np.pi)/dphi))
                    rp_withcatching[j] += np.fabs(qgphi_del[0,i])*qgphi_del[1,i]

        rp[0] += rp[n-1]
        rp[n-1] = rp[0]
        if polar:
            phi[n-1] = phi[0]

        rpint = 0
        for i in np.arange(n):
            rpint += rp[i]
            if rp_withcatching is not None:
                rpint += rp_withcatching[i]
        if rpint!=0:
            for i in np.arange(n):
                rp[i] = rp[i]*2*np.pi/(dphi*rpint)
                if rp_withcatching is not None:
                    rp_withcatching[i] = rp_withcatching[i]*2*np.pi/(dphi*rpint)
        
        if not polar:
            x_ticks = [-np.pi, -0.5*np.pi, 0, 0.5*np.pi, np.pi]
            labels = [r'$-\pi$', r'$-\pi/2}{2}$',r'$0$',r'$\pi/2$',r'$\pi$']
            plt.xticks(x_ticks, labels)
        plt.plot(phi,rp,color=c[ci])
        if rp_withcatching is not None:
            plt.plot(phi,rp_withcatching,'--',color=c[ci])
        ci+=1

    if save2 is not None:
        plt.savefig(save2)


def spectrum(t=None, particles='geip', colors='bgmrcyk', sptype='simple', axis=[], save2=None, data_folder=None, smooth=False, xlim=None, ylim=None,
        smooth_start=20, smooth_max=1000, smooth_width=50, window_type='triangular', multi_mev_threshold=None, clf=False, **kwargs):
    'spectrum() # plots spectrum for all particles\n\
    at t_end\n\
    Examples:\n\
    spectrum(10,sptype=\'loglog\') # energy distribution\n\
    in log-log axes.\n\
    spectrum(5,sptype=\'lin\',axis=[0,1,-2,2]) \n\
    # linear scaling of x and y axes with x from 0 to 1 and y from -2 to 2.\n'

    def window(i, width, type='square'):
        'Window function with specified width'
        if type=='square':
            if i>width or i<-width:
                return 0.0
            else:
                return 1.0 / (2*width + 1)
        elif type=='triangular':
            if i>width or i<-width:
                return 0.0
            elif i>=0:
                return (1.0 - 1.0*i/(width+1)) / (width+1)
            else:
                return (1.0 + 1.0*i/(width+1)) / (width+1)

    def smooth_array(x, smooth_start, smooth_max, max_width, wtype):
        'Returnes smoothed array using specified window, starting from smooth_start.\n\
        Window width linearly increases until x[i]==smooth_max and after that remains const=max_width'
        result = np.zeros(len(x))
        for i in np.arange(len(x)):
            width = 0
            
            if i > smooth_start and i < smooth_max:
                width = int(max_width * (i - smooth_start) * 1.0 / (smooth_max - smooth_start))
            elif i > smooth_start:
                width = max_width

            for j in np.arange(-width, width+1):
                if i+j>=0 and i+j<len(x):
                    result[i] += x[i+j] * window(j, width, wtype) 
        return result
    
    if clf:
        plt.clf()
    plt.xlabel('kinetic energy, MeV')
    if sptype=='energy' or sptype=='loglog':
        plt.ylabel(r'$\varepsilon dN/d\varepsilon$, a.u.')
    else:
        plt.ylabel(r'$dN/d\varepsilon$, a.u.')
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if axis!=[]:
        plt.axis(axis)
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    if t==None:
        t = resread.t_end - resread.dt
    resread.t = '%g' % t
    s = []
    deps = []
    ci = []
    i = -1
    for p in list(particles):
        i+=1
        if p=='e':
            s.append('')
            deps.append(resread.deps)
            ci.append(i%len(colors))
        elif p=='p':
            s.append('_p')
            deps.append(resread.deps_p)
            ci.append(i%len(colors))
        elif p=='g':
            s.append('_ph')
            deps.append(resread.deps_ph)
            ci.append(i%len(colors))
        elif p=='i':
            for cmr in resread.icmr:
                s.append('_'+str(cmr)+'_')
                deps.append(resread.deps_i)
                ci.append(i%len(colors))
    ret_val = {}
    for i in np.arange(len(s)):
        sp = resread.t_data('spectrum'+s[i]+resread.t,deps[i])
        if (resread.catching or (resread.dump_photons and s[i] == '_ph')) and include_deleted:
            for t1 in np.arange(0, float(resread.t)+0.001, resread.output_period):
                sp += resread.t_data('spectrum_deleted{0}{1:g}'.format(s[i], t1), deps[i])
        if multi_mev_threshold is not None:
            ret_val[s[i]] = {}
            multi_mev = 0.0
            for j in np.arange(len(sp[:,0])):
                if sp[j,0] > multi_mev_threshold:
                    multi_mev = multi_mev + sp[j,1] * resread.deps
                    
            print ('Number of multi-MeV particles of type [{0}] (W > {1} MeV) = {2}'.format(s[i], multi_mev_threshold, multi_mev)) 
            ret_val[s[i]]['multi_mev'] = multi_mev
            if s[i] == '' or s[i] == '_p':
                print ('Total charge = {0} nC'.format(multi_mev * 1.6e-10))
                ret_val[s[i]]['total_charge'] = multi_mev * 1.6e-10

        if smooth:
            sp[:,1] = smooth_array(sp[:,1], smooth_start, smooth_max, smooth_width, window_type)

        if sptype != 'lin':
            plt.yscale('log')
        if sptype=='energy' or sptype=='loglog':
            for j in np.arange(len(sp[:,0])):
                sp[j,1] = sp[j,0]*sp[j,1]
        if sptype=='loglog':
            plt.xscale('log')
        plt.plot(sp[:,0],sp[:,1],colors[ci[i]])

    if save2 is not None:
        plt.savefig(save2)
    if multi_mev_threshold is not None:
        return ret_val


directivity = 0
directivity_lat = 0
directivity_lng = 0
directivity_lat1 = 0
directivity_lng1 = 0
directivity_lat2 = 0
directivity_lng2 = 0


def mollweide(t=None, nlongitude=80, nlatitude=40, Nlevels=15, save2=None, data_folder=None, data=None, clf=False, **kwargs):
    '''Plots photon radiation pattern in Mollweide projection and computes (antenna-like) directivity.
    
    In the Mollweide projection, z-axis sticks out of the north pole and y-axis sticks out of the\n\
    projection center.

    Return value:
        lng, lat, rp - grid points and computed radiation pattern
    Arguments:
        t -- the time at which the radiation pattern is calculated
        data_folder (or df) -- folder to take the data from
        nlongitude, nlatitude -- determines the plot resolution
        data -- if specified, the other parameters are ignored and 'data' is plotted.
            It can be in one of the 2 formats: either (lng, lat, rp) as returned by this function, or just rp (2D array)
    Note:
        The function inits the global variables 'directivity', 'directivity_lat', etc.
    '''
    global directivity, directivity_lat, directivity_lng,\
    directivity_lat1, directivity_lng1, directivity_lat2,\
    directivity_lng2

    if data is None:
        data_folder = __get_data_folder(data_folder, kwargs)
        if data_folder is not None:
            resread.data_folder = data_folder
        resread.read_parameters()
        if t==None:
            t = resread.t_end - resread.dt
        resread.t = '%g' % t
        #
        lat = np.linspace(-np.pi/2,np.pi/2,nlatitude)
        lng = np.linspace(-np.pi,np.pi,nlongitude)
        rp = np.zeros((nlatitude,nlongitude))
        #
        qgthphi = resread.particles('phasespace_ph',['q','g','theta','phi'])
        print ('Number of particles: ' + str(len(qgthphi[0,:])))
        cnt = 0
        print('Computing radiation pattern...', end='')
        for i in np.arange(len(qgthphi[0,:])):
            if i % (len(qgthphi[0,:]) // 9 + 1) == 0:
                cnt+=1
                print ('{0}%..'.format(cnt*10), end='')
            # theta - latitude, from xy-plane, [-pi/2,pi/2]
            # phi - longitude, in xy-plane, (-pi,pi]
            y1 = (qgthphi[2,i]+np.pi/2)*(nlatitude-1)/np.pi
            j = np.floor(y1)
            y1 = y1 - j
            y2 = 1 - y1
            x1 = (qgthphi[3,i]+np.pi)*(nlongitude-1)/(2*np.pi)
            k = np.ceil(x1)
            x1 = k - x1
            x2 = 1 - x1
            if j!=nlatitude:
                rp[j,k] += y1*x1*qgthphi[0,i]*qgthphi[1,i]
                rp[j+1,k] += y2*x1*qgthphi[0,i]*qgthphi[1,i]
                rp[j+1,k-1] += y2*x2*qgthphi[0,i]*qgthphi[1,i]
                rp[j,k-1] += y1*x2*qgthphi[0,i]*qgthphi[1,i]
            else:
                rp[j,k] += y1*x1*qgthphi[0,i]*qgthphi[1,i]
                rp[j,k-1] += y1*x2*qgthphi[0,i]*qgthphi[1,i]
        print('100%')
        tmp1 = 0
        tmp2 = 0
        for i in np.arange(nlongitude):
            tmp1 += rp[0,i]
            tmp2 += rp[nlatitude-1,i]
        for i in  np.arange(nlongitude):
            rp[0,i] = tmp1
            rp[nlatitude-1,i] = tmp2
        for i in np.arange(1,nlatitude-1):
            rp[i,0] += rp[i,nlongitude-1]
            rp[i,nlongitude-1] = rp[i,0]
        for j in np.arange(nlongitude):
            for i in np.arange(1,nlatitude-1):
                rp[i,j] = rp[i,j]/np.cos(lat[i])
            rp[0,j] = rp[0,j]*8/np.cos(lat[1])/(nlongitude-1)
            rp[nlatitude-1,j] = rp[nlatitude-1,j]*8/np.cos(lat[1])/(nlongitude-1)
        integral = 0
        integral1 = 0
        for j in np.arange(nlongitude):
            for i in np.arange(1,nlatitude-1):
                integral += rp[i,j]*np.cos(lat[i])
                integral1 += np.cos(lat[i])
        integral += rp[0,0]*np.cos(lat[1])/8
        integral += rp[nlatitude-1,0]*np.cos(lat[1])/8
        integral1 += np.cos(lat[1])/8
        integral1 += np.cos(lat[1])/8
        for j in np.arange(nlongitude):
            for i in np.arange(nlatitude):
                rp[i,j] = rp[i,j]*integral1/integral
        directivity = rp.max()
        directivity_lat = lat[np.floor(rp.argmax()/nlongitude)]
        directivity_lng = lng[rp.argmax()%nlongitude]
        directivity_lat1 = lat[np.floor(rp[:,:nlongitude/2].argmax()/(nlongitude/2))]
        directivity_lng1 = lng[rp[:,:nlongitude/2].argmax()%(nlongitude/2)]
        directivity_lat2 = lat[np.floor(rp[:,nlongitude/2:].argmax()/(nlongitude-nlongitude/2))]
        directivity_lng2 = lng[nlongitude/2 + rp[:,nlongitude/2:].argmax()%(nlongitude/2)]
    else:
        if len(data) == 3: # lng, lat, rp
            lng, lat, rp = data
        else:              # only rp
            rp = data
            lat = np.linspace(-np.pi/2,np.pi/2,data.shape[0])
            lng = np.linspace(-np.pi,np.pi,data.shape[1])
    if clf:
        plt.clf()
    f = plt.figure()
    ax = f.add_subplot(111,projection='mollweide')
    ax.contour(lng,lat,rp,Nlevels,origin=None)

    if save2 is not None:
        plt.savefig(save2)
    return lng, lat, rp


def field(t=0, field='ex', plane='xy', field2=None, fmax=None, data_folder=None, extent=None, xlim=None, ylim=None, axis=[], save2=None, clf=False, **kwargs):
    '''Plots fields in the specified plane.

    Arguments:
        t -- the time at which the field is plotted
        field -- the field to be plotted
        field2 (optional) -- the second field which can be added to the 1st one.
            Example: field='ey', field2='-bz' plots Ey-Bz value
        plane -- either 'xy', 'xz' or 'yz'. The 3rd coordinate is pre-determined (set in the Quill config file)
        data_folder (or df) -- folder to take data from
    '''
    resread.t = '%g' % t
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    #
    if clf:
        plt.clf()
    title = tex_format(field)
    if field2 is not None:
        if field2[0] == '-':
            title += '-' + tex_format(field2[1:])
        else:
            title += '+' + tex_format(field2)
    plt.title(title)
    plt.xlabel(tex_format(plane[0]) + '$/\lambda$')
    plt.ylabel(tex_format(plane[1]) + '$/\lambda$')
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if plane=='xz':
        xlength = resread.nx*resread.dx
        ylength = resread.nz*resread.dz
    elif plane=='yz':
        xlength = resread.ny*resread.dy
        ylength = resread.nz*resread.dz
    else:
        xlength = resread.nx*resread.dx
        ylength = resread.ny*resread.dy
    if axis!=[]:
        plt.axis(axis)
    #
    f = resread.density(field,plane)

    # additive second field
    if field2 is not None:
        if field2[0] == '-':
            f2 = -1.0 * resread.density(field2[1:],plane)
        else:
            f2 = resread.density(field2,plane)
        f += f2

    if fmax==None:
        fmax = np.max( [-np.min(f), np.max(f)] )
        print('fmax = {0}'.format(fmax))
    #
    if extent == None:
        extent = (0,xlength,0,ylength)
    elif extent == 'xi':
        extent = ( resread.xi( 0, t ), resread.xi( xlength, t ), 0, ylength )
    plt.imshow(f,'bwr',interpolation='none',vmin=-fmax,vmax=fmax,origin='lower',extent=extent)
    #
    if save2 is not None:
        plt.savefig(save2)


def energy(data_folder=None, save2=None, include_deleted=True, clf=False, species='epgiwt', **kwargs):
    'Plots energy of electrons, ions, etc. vs time'
    'species - what should be plotted: e, p, g, i - particles, w - field energy, t - total'
    #Important: several types of ions are not supported

    def safe_sum(a, b):
        result = np.array([])
        min_len = min(len(a), len(b))
        for i in np.arange(min_len):
            result = np.append(result, a[i] + b[i])
        return result, min_len

    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    tmp = resread.t_data('energy')
    if clf:
        plt.clf()
    if 'w' in species:
        plt.plot(tmp[:,0], tmp[:,1], 'k') # em fields
    if 'e' in species:
        plt.plot(tmp[:,0], tmp[:,2], 'g') # electrons
    if 'p' in species:
        plt.plot(tmp[:,0], tmp[:,3], 'r') # positrons
    if 'g' in species:
        plt.plot(tmp[:,0], tmp[:,4], 'b') # hard photons
    total_energy = np.sum(tmp[:,1:5], axis=1)
    if (resread.n_ion_populations>0):
        total_energy += tmp[:,5]
        if 'i' in species:
            plt.plot(tmp[:,0], tmp[:,5], 'm') # ions
    if 't' in species:
        plt.plot(tmp[:,0], total_energy,'--k')

    if (resread.catching or (resread.dump_photons and 'g' in species)) and include_deleted:
        # deleted energy
        tmp_del = resread.t_data('energy_deleted')
        if 'e' in species:
            sum_e, min_len = safe_sum(tmp_del[:,1], tmp[:,2])
            plt.plot(tmp_del[:min_len,0], sum_e, '--g') # electrons
        if 'p' in species:
            sum_p, min_len = safe_sum(tmp_del[:,2], tmp[:,3])
            plt.plot(tmp_del[:min_len,0], sum_p, '--r') # positrons
        if 'g' in species:
            sum_g, min_len = safe_sum(tmp_del[:,3], tmp[:,4])
            plt.plot(tmp_del[:min_len,0], sum_g, '--b') # hard photons
        
        total_del_energy, min_len = safe_sum(total_energy, np.sum(tmp_del[:,1:4], axis=1))
        if (resread.n_ion_populations>0):
            total_del_energy, min_len = safe_sum(total_del_energy, tmp_del[:,4])
            sum_i, min_len = safe_sum(tmp_del[:,4], tmp[:,5])
            if 'i' in species:
                plt.plot(tmp_del[:min_len,0], sum_i, '--m') # ions
        if 't' in species:
            plt.plot(tmp[:min_len,0], total_del_energy, ':k') # sum energy
    
    plt.xlabel('$ct/\lambda$')
    plt.ylabel('Energy, J')

    if save2 is not None:
        plt.savefig(save2)


def tracks2(space=['x', 'y'], tracks=None, particles='e', save2=None, data_folder=None, filter=None, clf=True, **kwargs):
    'Plots 2d tracks in the specified *space*'
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    if not tracks:
        if len(particles) > 1 or not particles[0] in 'epg':
            raise ValueError('tracks2() does not support multiple particle species, as well as ions')
        tracks = resread.tracks(particles, filter)
    if clf:
        plt.clf()
    if len(space) == 2:
        x = space[0]
        y = space[1]
        plt.xlabel(tex_format(x))
        plt.ylabel(tex_format(y))
        for track in tracks:
            plt.plot(track[x], track[y])
    elif len(space) == 3:
        x = space[0]
        y = space[1]
        c = space[2]
        plt.xlabel(tex_format(x))
        plt.ylabel(tex_format(y))
        plt.title(tex_format(c))
        for track in tracks:
            plt.scatter(track[x], track[y], c = track[c], s = 1, edgecolor='none')
        plt.colorbar()
    if save2 is not None:
        plt.savefig(save2)


def N(data_folder=None, particles='gep', save2=None, clf=False, **kwargs):
    'Plots number of particles over time. Ions are not currently supported'
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters();
    a = resread.t_data('N')
    t = a[:,0]
    Ne = a[:,1]
    Np = a[:,2]
    Ng = a[:,3]

    lw = 1.5 # linewidth
    if clf:
        plt.clf()
    if 'e' in particles:
        plt.plot(t, Ne, 'g--', linewidth = lw, label = r'$N_e$')
    if 'p' in particles:
        plt.plot(t, Np, 'r-', linewidth = lw, label = r'$N_p$')
    if 'g' in particles:
        plt.plot(t, Ng, 'b:', linewidth = lw, label = r'$N_g$')

    plt.legend(loc = 'best', fontsize = 'medium')
    plt.xlabel(r'$ct/\lambda$')
    plt.ylabel(r'$N$')
    plt.yscale('log')

    if save2 is not None:
        plt.savefig(save2)


def onaxis(t, particles='we', colors='rgbcmyk', norm='true',
        data_folder=None, save2=None, plotargs={}, rrargs={}, clf=False, **kwargs):
    'Plot of particle density, fields etc. along the x-axis.\n\
    *norm* can be set to \'optimal\', \'true\' or to array of desired maximal values.\n\
    For rrargs see help(qplot.resread.onaxis),\n\
    for plotargs see help(plt.plot).\n\
    \n\
    Examples:\n\
    qplot.onaxis(0, \'e\'),\n\
    qplot.onaxis(0, \'ep\', norm = \'true\'),\n\
    qplot.onaxis(0, \'we\', norm = [0.5, 1]),\n\
    qplot.onaxis(15,[\'ey\', \'bz\', \'e\'], \'rgb\', plotargs = {linewidth: 1.2})\n\
    qplot.onaxis(4, \'g\', rrargs = {\'av\': \'y\'}),\n\
    qplot.onaxis(4, \'p\', rrargs = {\'sx\': 10, \'sz\': 3}).'
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    x = resread.onaxis('x', **rrargs)
    a = []
    for p in particles:
        if p == 'e':
            filename = 'rho'
        elif p == 'g':
            filename = 'rho_ph'
        elif p == 'p':
            filename = 'rho_p'
        elif p == 'i':
            filename = 'irho_' + str(resread.icmr[0]) + '_'
        else:
            filename = p
        resread.t = '%g' % t
        if p == 'e':
            a.append(-resread.onaxis(filename, **rrargs))
        else:
            a.append(resread.onaxis(filename, **rrargs))
    ma = 0
    for i, b in enumerate(a):
        mb = max(b)
        print('max value for', particles[i], '=', mb)
        if mb > ma:
            ma = mb
    if norm == 'optimal':
        for i, b in enumerate(a):
            mb = max(b)
            if mb < ma / 4:
                a[i] = a[i] / mb * ma / 4
    elif norm != 'true':
        for i, m in enumerate(norm):
            a[i] = norm[i] * a[i] / max(a[i])

    if clf:
        plt.clf()
    for i, b in enumerate(a):
        plt.plot(x, b, color = colors[i % len(colors)], **plotargs)
    plt.xlabel('$x$')
    if save2 is not None:
        plt.savefig(save2)


def ion():
    """
    Invokes 'plt.ion()'
    """
    plt.ion()


def ioff():
    """
    Invokes 'plt.ioff()'
    """
    plt.ioff()


def reset_style():
    mpl.rcParams.update(rc_backup) 

import __main__
if not plt.isinteractive() and not hasattr(__main__, '__file__'):
    print("For the interactive regime use IPython in the Pylab mode ('ipython --pylab' or '%pylab'), 'plt.ion()', "
          "or 'qplot.ion()'")

lw = 0.7 # linewidth for border
lwl = 1.0 # linewidth for lines in plots
font = {'family' : 'serif', 'serif' : 'cmr10', 'size' : 9}
rc_backup = mpl.rcParams.copy()
mpl.rc('font', **font)
mpl.rc('lines', linewidth=lwl)
mpl.rc('axes', linewidth=lw)
mpl.rc('axes', unicode_minus=False)
mpl.rc('figure', figsize=(3.5,2.5), dpi=200, autolayout=True)
mpl.rc('mathtext' ,fontset='cm')
mpl.rc('savefig',dpi=300)

tcmap.red()
tcmap.green()
tcmap.blue()
tcmap.orange()
tcmap.purple()
