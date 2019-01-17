import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot
import matplotlib.animation as animation
import numpy as np
import re

rx = re.compile(r'[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?', re.VERBOSE)

def read_area(file_name, var):
    data = open(file_name, 'r')
    data.readline()
    t = data.readline().split('\t')[1][:-1]
    s = data.readline().split('\t')[1][:-1]
    [x1, y1, x2, y2, nx, ny] = [float(i) for i in rx.findall(data.readline())]
    nx, ny = int(nx), int(ny)
    data.readline()
    rp = []
    ip = []
    for i in data:
        rp.append(float(rx.findall(i)[var * 2 + 3]))
        ip.append(float(rx.findall(i)[var * 2 + 4]))
    return np.asarray(rp).reshape(nx, ny), np.asarray(ip).reshape(nx, ny)

def plot_area(file_name, var, phase = 0):
    rp, ip = read_area(file_name, var)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    w = rp * np.cos(phase) + ip * np.sin(phase)
    return plt.imshow(w,interpolation="none", cmap='Blues')

def read_circle(file_name, var):
    data = open(file_name, 'r')
    data.readline()
    t = data.readline().split('\t')[1][:-1]
    s = data.readline().split('\t')[1][:-1]
    [x, y, r, n] = [float(i) for i in rx.findall(data.readline())]
    n = int(n)
    data.readline()
    a = []
    rp = []
    ip = []
    for i in data:
        a.append(float(rx.findall(i)[2]))
        rp.append(float(rx.findall(i)[var * 2 + 3]))
        ip.append(float(rx.findall(i)[var * 2 + 4]))
    return np.asarray(a), np.asarray(rp), np.asarray(ip)

def anim_area(file_name, var):
    d = 10 * 2*np.pi / 400
    rp, ip = read_area(file_name, var)
    f = plt.figure()
    ax = plt.Axes(f, [0., 0., 1., 1.])
    ax.set_axis_off()
    f.add_axes(ax)
    im = plt.imshow(rp, interpolation='none', cmap='Blues', origin='lower')
    a = lambda i : (im.set_array(rp * np.cos(d*i) + ip*np.sin(d*i)),)
    return animation.FuncAnimation(f, a, interval=50, frames=400, blit=False)
