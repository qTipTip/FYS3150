import os
from numpy import pi
def run_gauss_legendre():
    n_values = range(5, 45, 5)
    limit_values = [3]

    exact = 5 * pi**2 / 16.0**2

    for n in n_values:
        for r in limit_values:
            os.system('./a.out %d 1 %d %d' % (n, -r, r))

def run_gauss_laguerre():
    n_values = range(5, 45, 5)
    limit_values = [3]

    exact = 5 * pi**2 / 16.0**2

    for n in n_values:
        for r in limit_values:
            os.system('./a.out %d 2 %d %d' % (n, -r, r))

def run_bruteforce_MC():
    n_values = range(10)
    for n in n_values:
        os.system('./a.out %d 3 -3 3' % 10**n)

def read_gauss_legendre():
    times = {}
    values = {}
     
    with open('../data/gauss_legendre.dat', 'r') as infile:
        for line in infile:
            line = line.split()
            n = int(line[0])
            t = float(line[2])
            result = float(line[1])
            a = int(line[3])
            b = int(line[4])

            if n in times:
                times[n] += t
            else:
                times[n] = t

            if n in values:
                values[n].append(result)
            else:
                values[n] = [result]

        # Averaging time elapsed and changing to seconds
        for key in times:
            times[key] /= 5 * 1.0e9

    return times, values

def dict_to_tex(dictionary, filename):
    output = '\\begin{tabular}{cc}\n'
    output += '$n$ & average time elapsed [s]\\\\\n'
    output += '\\hline\n'

    for key in dictionary:
        output += '%d & %0.4f \\\\\n' % (key, dictionary[key])
    output+= '\\end{tabular}\n'

    with open(filename, 'w') as outfile:
        outfile.write(output)

def matrix_to_tex(dictionary, filename):
    exact = 5 * pi * pi / 16.0**2
    output = '\\begin{tabular}{cccccc}\n'
    for key in dictionary:
        if key == 30:
            output += str(key) + "\\\\\n"
        else:
            output += str(key) + " & " 
    output += '\\hline\n'
    for key in dictionary:
        for val in dictionary[key][:-1]:
            output +=  "%0.4f & " % val - exact
        output += str(dictionary[key][-1])
        output += "\\\\\n"
    output += '\\end{tabular}\n'

    with open(filename, 'w') as outfile:
        outfile.write(output)
    
if __name__ == '__main__':
    #run_gauss_legendre()
    #run_gauss_laguerre()
    run_bruteforce_MC()
    # times, values = read_gauss_legendre()
    # dict_to_tex(times, 'times_legendre.tex') Done
    # matrix_to_tex(values, '../article_revised/values_legendre.tex')

