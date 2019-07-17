from transfer_matrix_method import *

wavelength = []
index = []
extinct = []

params = '/Users/garrek/projects/pistachio/data/Au_SiO2_n_Ciesielski.csv'

def get_data_from_csv(path):

    with open(path, 'r') as params:
        reader = csv.reader(params)
        next(reader, None)
        for row in reader:
            wl = float(row[0])
            n = float(row[1])
            wavelength.append(wl)
            index.append(n)
            if row[2]:
                K = float(row[2])
                extinct.append(K)


get_data_from_csv(params)
    
f = sp.interpolate.interp1d(wavelength, index)

num_points = 1000
n_min = 0.2
n_max = 15

new_x = np.linspace(n_min, n_max, num=num_points, endpoint=True)

plt.plot(wavelength, index, 'o', new_x, f(new_x), '-')
plt.xlim(0, 20)
plt.show()