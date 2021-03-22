import stochpy
import numpy as np
import mpmath
import pandas as pd
import matplotlib.pyplot as mplt

from bioservices import UniProt         # Import bioservices module, to run remote UniProt queries



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


def open_excel_file():
    i = 1
    gene_dict = {}
    kdp = 1.0 / 150
    file = pd.read_excel(r'F:\Studia\Gillespie_stoch_model\NIHMS211541-supplement-TableS6.xls')
    df = pd.DataFrame(file, columns=['Gene Name', 'Mean_RNAseq', 'Life time_RNAseq'])
    df = df.rename(columns={'Gene Name': 'Gene', 'Life time_RNAseq':'Lifetime'})
    for index, rows in df.iterrows():
        if rows.isnull().values.any() == False:
            gene_dict[rows.Gene] = [rows.Mean_RNAseq, rows.Lifetime]
    return gene_dict


def get_sequences_lengths(gene_dict):           # Get length of protein (E.Coli) from first found hit in UniProt database.
    service = UniProt()
    length = []
    for gene in gene_dict.keys():
        query = gene                            # Build a query string
        result = service.search(query)          # Send the query to UniProt, and catch the search result in a variable
        result = result.split("\n")
        for i in result:
            i = i.split("\t")
            tmp = i[5].split()
            if "coli" in tmp:
                length.append(i[6] + ",")
                print i[6]
                break
    return length


def count_parameters(translation_rate):       # Do poprawy
    kp = float(rows.lengthaa) / (kp_av_in_s * 60)
    kdm = 1 / float(rows.Lifetime)
    km = float(rows.Mean_RNAseq) * kdm
    gene_dict[rows.Gene] = [km, kp, kdm, kdp]


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
    stochpy.plt.savefig("F:\Studia\Gillespie_stoch_model\PlotSpecies for gene\PlotSpecies %s.png" % gene_name, format='png')
    stochpy.plt.show()


    smod2.PlotSpeciesDistributions(species2plot='P')
    stochpy.plt.step(x, y, color='r')
    stochpy.plt.title("Histogram of protein, gene %s" % gene_name)
    stochpy.plt.savefig("F:\Studia\Gillespie_stoch_model\Histogram of protein for genes\Histogram %s.png" % gene_name, format='png')
    stochpy.plt.show()


def ergodicity_check(gene_name, data_set):                          #Do poprawy
    smod = stochpy.SSA(IsInteractive=False)
    smod.Model("Double_step.psc")
    smod.Model(model_file="Double_step.psc", dir="C:\Stochpy\pscmodels")
    smod.DoStochSim(method="direct", trajectories=3000, mode="time", end=100)

    trajectories_number = smod.data_stochsim.simulation_trajectory
    smod.PlotSpeciesTimeSeries()
    stochpy.plt.show()
    stochpy.plt.savefig("Wykres.png")
    smod.ShowOverview()

    last_val = Last_values_ofTrajectories(trajectories_number, smod)  # Proba plotowania histogramu po osttanich wartosciach w trajektoriach.
    x, y = percent(last_val)
    smod.PlotSpeciesDistributions(species2plot='P')
    stochpy.plt.step(x,y, color='r')
    stochpy.plt.show()



data = open_excel_file()
print len(data.keys())
length = get_sequences_lengths(data)
file = open('F:\Studia\Gillespie_stoch_model\lengths.txt', 'w')
file.writelines(length)

"""
for gene in data.keys():
    print(gene, data[gene])
    simple_simulation(gene, data)
"""
import io

# Import Seaborn for graphics and plotting
import seaborn as sns



