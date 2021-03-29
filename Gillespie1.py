import stochpy
import numpy as np
import mpmath
import pandas as pd
import matplotlib.pyplot as mplt
from bioservices import UniProt         # Import bioservices module, to run remote UniProt queries


def GammaDis(km, kp, kdm, kdp, max_val):         #Funkcja do obliczania rozkladu Gamma
    xv = np.linspace(0,max_val,1000)
    yv = []
    a = km/kdp
    b = kp/kdm
    for x in xv:
        p = x**(a-1)*mpmath.exp(-x/b)/((b**a)*mpmath.mp.gamma(a))
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
    gene_dict = {}
    gene_list = {}
    file = pd.read_excel(r'F:\Studia\Gillespie_stoch_model\NIHMS211541-supplement-TableS6.xls')
    df = pd.DataFrame(file, columns=['Gene Name', 'Mean_RNAseq', 'Life time_RNAseq', 'Mean_Protein'])
    df = df.rename(columns={'Gene Name': 'Gene', 'Life time_RNAseq':'Lifetime'})
    for index, rows in df.iterrows():
        if rows.isnull().values.any() == False:
            gene_dict[rows.Gene] = [rows.Mean_RNAseq, rows.Lifetime, rows.Mean_Protein]
            gene_list[rows.Gene] = True
        else:
            gene_list[rows.Gene] = False
    return gene_dict, gene_list


def open_sequences_length_file():               # Opening existing file which contains protein sequences lengths from UniProt.
    f = open('F:\Studia\Gillespie_stoch_model\lengths.txt', 'r')
    lengths = f.read().split(",")[:-1]
    return lengths


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
                break
    file = open('F:\Studia\Gillespie_stoch_model\lengths.txt', 'w')
    file.writelines(length)
    return length


def count_parameters(gene_dict, sequences_lengths, average_translation_rate, kp1 = False, kp2 = False, kp3 = False):           #Funkcja pomocnicza pomagajaca policzyc parametry symulacji km,kp,kdm,kdp.
    parameters_dict = {}
    gene_length_id = 0
    for gene in gene_dict.keys():
        if kp1 is True:
            kp = float(sequences_lengths[gene_length_id]) / (average_translation_rate * 60)
            print kp
        elif kp2 is True:
            kp = 60.0
        else:
            kp = (float(gene_dict[gene][2])/float(gene_dict[gene][0]))*float(gene_dict[gene][1])
            print("%", kp)
        kdm = 1 / float(gene_dict[gene][1])
        km = float(gene_dict[gene][0]) * kdm
        kdp = 1.0/150
        parameters_dict[gene] = [km, kp, kdm, kdp]
        gene_length_id += 1
    return parameters_dict


def kp_to_excel(gene_parameters_list, genel_list, kp):
    column = []
    i = 0
    file = pd.read_excel(r'F:\Studia\Gillespie_stoch_model\NIHMS211541-supplement-TableS6.xls')
    df = pd.DataFrame(file)
    for val in genel_list.keys():
        if genel_list[val] is False:
            column.append("")
        else:
            column.append(gene_parameters_list[val][1])
            i += 1
    df[kp] = column
    df.to_excel('F:\Studia\Gillespie_stoch_model\NIHMS211541-supplement-TableS6.xls')


def simple_simulation(gene_name, data_set):                         #Funkcja przeprowadzajaca 1 symulacje dla danego genu.
    smod2 = stochpy.SSA(IsInteractive=False)
    smod2.Model("Double_step.psc")
    smod2.Model(model_file="Double_step.psc", dir="C:\Stochpy\pscmodels")
    smod2.ChangeParameter("km", data_set[gene_name][0])
    smod2.ChangeParameter("kp", data_set[gene_name][1])
    smod2.ChangeParameter("kdm", data_set[gene_name][2])
    smod2.ChangeParameter("kdp", data_set[gene_name][3])
    burst_size = float(data_set[gene_name][1])/float(data_set[gene_name][2])

    if burst_size >= 6:
        smod2.DoStochSim(method="direct", trajectories=1, mode="time", end=18000)
    else:
        smod2.DoStochSim(method="direct", trajectories=1, mode="time", end=9000)

    max_val = 0
    for i in smod2.data_stochsim.getAllSimData():
        if i[3] >= max_val:
            max_val = i[3]
    max_val = int(mpmath.ceil(max_val))

    x, y = GammaDis(data_set[gene_name][0], data_set[gene_name][1], data_set[gene_name][2], data_set[gene_name][3], max_val)
    smod2.PlotSpeciesTimeSeries()
    stochpy.plt.title("PlotSpecies, gene %s, b=%s" % (gene_name, burst_size))
    stochpy.plt.savefig("F:\Studia\Gillespie_stoch_model\PlotSpecies for genes, kp2\PlotSpecies %s.png" % gene_name, format='png')
    stochpy.plt.show()

    smod2.PlotSpeciesDistributions(species2plot='P')
    stochpy.plt.step(x, y, color='r')
    stochpy.plt.title("Histogram of protein, gene %s, b=%s" % (gene_name, burst_size))
    stochpy.plt.savefig("F:\Studia\Gillespie_stoch_model\Histogram of protein, kp2\Histogram %s.png" % gene_name, format='png')
    stochpy.plt.show()


def ergodicity_check(gene_name, data_set):                      #Do poprawy
    smod = stochpy.SSA(IsInteractive=False)
    smod.Model("Double_step.psc")
    smod.Model(model_file="Double_step.psc", dir="C:\Stochpy\pscmodels")
    smod.ChangeParameter("km", data_set[gene_name][0])
    smod.ChangeParameter("kp", data_set[gene_name][1])
    smod.ChangeParameter("kdm", data_set[gene_name][2])
    smod.ChangeParameter("kdp", data_set[gene_name][3])
    smod.DoStochSim(method="direct", trajectories=3000, mode="time", end=100)

    smod2 = stochpy.SSA(IsInteractive=False)
    smod2.Model("Double_step.psc")
    smod2.Model(model_file="Double_step.psc", dir="C:\Stochpy\pscmodels")
    smod2.ChangeParameter("km", data_set[gene_name][0])
    smod2.ChangeParameter("kp", data_set[gene_name][1])
    smod2.ChangeParameter("kdm", data_set[gene_name][2])
    smod2.ChangeParameter("kdp", data_set[gene_name][3])
    smod2.DoStochSim(method="direct", trajectories=1, mode="time", end=3000)

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


average_translation_rate = 8.0
data, genel_list = open_excel_file()     # Read data from excel table

sequences_length = open_sequences_length_file()                 # Read data from sequence's lengths file
gene_parameters_dict1 = count_parameters(data, sequences_length, average_translation_rate, True, False, False)
gene_parameters_dict3 = count_parameters(data, sequences_length, average_translation_rate, False, False, True)       # Create dictionary that contains name of gene and parameters for simulation


kp_to_excel(gene_parameters_dict1, genel_list, 'kp1')                 # Export values of kp1 to S6 Table
kp_to_excel(gene_parameters_dict3, genel_list, 'kp3')
kp1 = []
kp3 = []


for gene in gene_parameters_dict1.keys():
    kp1.append(gene_parameters_dict1[gene][1])
for gene in gene_parameters_dict3.keys():
    kp3.append(gene_parameters_dict3[gene][1])
kp2 = [60.0]*len(kp1)
print(len(kp2))
mplt.scatter(kp1, kp3)
mplt.plot(kp1, kp2, color='r')
mplt.savefig("F:\Studia\Gillespie_stoch_model\Scatterplot k1 k3.png")
mplt.show()

#for gene in gene_parameters_dict.keys():                        # Do simulation for all genes
 #   simple_simulation(gene, gene_parameters_dict)



