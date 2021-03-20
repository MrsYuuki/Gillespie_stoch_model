import stochpy
import numpy as np
import mpmath
import matplotlib.pyplot as mplt
from matplotlib.ticker import PercentFormatter

def GammaDis(km, kp, kdm, kdp):         #Funkcja do obliczania rozkladu Gamma
    xv = np.linspace(0,400,1000)
    yv = []
    a = km/kdp
    b = kp/kdm
    for x in xv:
        p = x**(a-1)*np.exp(-x/b)/((b**a)*mpmath.mp.gamma(a))
        yv.append(p)
    print "Srednia czestosc burstow: ",a, " Sredni rozmiar burstow: ",b
    return xv, yv


def Last_values_ofTrajectories(trajectories, model):       #Funkkcja wydobywajace ostatnie wartosci natezenia bialek w trajektorii.
    last_trajectory_values = []
    for trajectory in range(1,trajectories+1):
        model.GetTrajectoryData(trajectory)
        last_trajectory_values.append(model.data_stochsim.getAllSimData()[-1][3])
    return last_trajectory_values


def percent(last_values):
    length = len(last_values)
    x_value = []
    y_value = []
    tmp_dict = {}
    for value in last_values:
        if value not in tmp_dict.keys():
            tmp_dict[value] = 1
        else:
            tmp_dict[value] += 1
    for values in tmp_dict:
        x_value.append(values)
        y_value.append(tmp_dict[values]/float(length))
    return x_value,y_value


def open_file(file_name, kp_av_in_s):
    file = open(file_name, "r")
    file_lines = {}
    line = file.readline()
    while line!="":
      n_line = line.replace("\n", "").split()
      kdp = 1.0/150
      kp = float(n_line[5])/(kp_av_in_s*60)
      kdm = 1/float(n_line[2])
      km = float(n_line[1]) * kdm
      file_lines[n_line[0]]=[km, kp, kdm, kdp]
      line = file.readline()
    return file_lines


def simple_simulation(gene_name, data_set):
    smod2 = stochpy.SSA(IsInteractive=False)
    smod2.Model("Double_step.psc")
    smod2.Model(model_file="Double_step.psc", dir="C:\Stochpy\pscmodels")
    smod2.ChangeParameter("km", data_set[gene_name][0])
    smod2.ChangeParameter("kp", data_set[gene_name][1])
    smod2.ChangeParameter("kdm", data_set[gene_name][2])
    smod2.ChangeParameter("kdp", data_set[gene_name][3])
    smod2.DoStochSim(method="direct", trajectories=1, mode="time", end=9000)
    x, y = GammaDis(data_set[gene_name][0], data_set[gene_name][1], data_set[gene_name][2], data_set[gene_name][3])

    smod2.PlotSpeciesTimeSeries()
    stochpy.plt.title("PlotSpecies, gene %s" % gene_name)
    stochpy.plt.show()

    smod2.PlotSpeciesDistributions(species2plot='P')
    stochpy.plt.step(x, y, color='r')
    stochpy.plt.title("Histogram of protein, gene %s" % gene_name)
    stochpy.plt.show()


data_set = open_file("C:\Users\Oliwia\PycharmProjects\pythonProject2\dane.txt", 8.4)
print(data_set)

simple_simulation("apaG", data_set)
"""smod = stochpy.SSA(IsInteractive=False)
smod.Model("Double_step.psc")
smod.Model(model_file="Double_step.psc", dir="C:\Stochpy\pscmodels")
smod.DoStochSim(method="direct", trajectories=3000, mode="time", end=100)



trajectories_number = smod.data_stochsim.simulation_trajectory
smod.PlotSpeciesTimeSeries()
stochpy.plt.show()
stochpy.plt.savefig("Wykres.png")
smod.ShowOverview()

last_val = Last_values_ofTrajectories(trajectories_number)         #Proba plotowania histogramu po osttanich wartosciach w trajektoriach.
x1,y1 = percent(last_val)

x,y = GammaDis(50,20,4,1)

#smod.PlotSpeciesDistributions(species2plot='P')                   #Dofitowanie wykresu rozkladu Gamma do histogramu ilosci bialka.
#stochpy.plt.step(x,y, color='r')
#stochpy.plt.show()"""







