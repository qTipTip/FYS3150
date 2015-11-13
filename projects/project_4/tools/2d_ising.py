import numpy as np
from PIL import Image
from random import randint, random
import shutil, glob, os

def create_image(configuration, dim_x, dim_y):
    """
    Given a matrix representing the current configuration of spins - assigns
    the color blue to a pixel with positive sign and color white to pixel with
    negative sign.
    """
    image = Image.new('RGB', (dim_x, dim_y), 'black')
    pixels = image.load()

    for a, i in enumerate(configuration):
        for b, j in enumerate(i):
            if j == 1:
                pixels[a, b] = (255, 0, 0)
            else:
                pixels[a, b] = (0,0, 255)

    return image

def initialize(dim_x, dim_y):
    """
    given the dimensions of a grid
    initializes it with random spin
    """

    matr = np.random.randint(2, size=(dim_x, dim_y))
    for a, i in enumerate(matr):
        for b, j in enumerate(i):
            if j == 0:
                matr[a, b] = -1
    return matr
def metropolis():
    i, j = [randint(0, dim_x-1), randint(0, dim_y-1)]
    E = get_energy(matr, i, j, dim_x, dim_y)
    if E > 0:
        matr[i, j] = -matr[i,j]
    else:
        r = random()
        if r < np.exp(2 * E / T):
            matr[i, j] = -matr[i,j]

def get_energy(configuration, i, j, dim_x, dim_y):

    this = configuration[i, j] 

    up =  configuration[i, j-1] if (j > 0) else configuration[i, dim_y-1]
    down =  configuration[i, j+1] if (j < dim_y-1) else configuration[i, 0]
    left = configuration[i-1, j] if (i > 0) else configuration[dim_x-1, j]
    right = configuration[i+1, j] if (i < dim_x-1) else configuration[0, j]

    return -this * (up + down + left + right)


if __name__ == '__main__':
    subdir = 'animation'
    if os.path.isdir(subdir): # if the directory is already present, remove it
        shutil.rmtree(subdir)
    os.mkdir(subdir) # create directory
    os.chdir(subdir) # enter new directory
    dim_x = 400
    dim_y = 400
    N = 30000000
    T_values = np.linspace(1, 5, N)
    matr = initialize(dim_x, dim_y)
    imag = create_image(matr, dim_x, dim_y)
    imag.save('start.png')
    
    framecounter = 0
    for T, n in zip(T_values, range(N)):
        metropolis()
        if n % 50000 == 0:
            print 'ping', n
            create_image(matr, dim_x, dim_y).save('tmp_%5d.png' % framecounter)
            framecounter += 1

    cmd = 'convert -delay 4 tmp_*.png movie.gif'
    os.system(cmd)

    for a in glob.glob('tmp*.png'):
        os.remove(a)
    os.chdir(os.pardir)
