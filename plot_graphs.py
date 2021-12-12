import os

ns = [10, 80, 1280, 20480]
ms = [  
        ("1", "#0533ed"),
        ("2", "#00FF00"),
        ("3A", "#FF0000"),
         ("3B", "#0000FF"),
        ("4A", "#777777"),
        ("4B", "#FF00FF"),
         #("5",  "#00FFFF")
    ]




for n in ns:
   
    s = f"""set terminal png size 1800,1350 enhanced font 'Verdana,9';
        set output 'n_f_plot_{n}.png';

        set size ratio 0.13;

        set format xy "%g";


        set grid;

        plot """



    filename = f'plot_{n}.plt'
    f = open(filename, 'w+')
    for (method, color) in ms:
        s += f''' "method_{method}_{n}.txt" title "{method}-{n}" with lines linecolor rgb '{color}' linetype 1 linewidth 1 ,'''

    f.write(s)
    f.close()
    os.system(f"gnuplot {filename}")

# filecontent = """set terminal png size 1800,1350 enhanced font 'Verdana,9';
# set output 'n_f_plot_10.png';

# set size ratio 0.13;

# set format xy "%g";


# set grid;

# plot  \
# "method_1_10.txt" title "1-5 function, f" with lines linecolor rgb '#0533ed' linetype 1 linewidth 1 ,\
#  "method_2_10.txt" title "2-5, g" with lines linecolor rgb '#00FF00' linetype 1 linewidth 1 ,\
#  #"points.txt" title "Points" with points linecolor rgb '#000000' #pointtype 6

# """






file = open("method_EPS_MAX.csv", "r")
v = []
for line in file:
    r = line.strip()
    r = r.strip(',')
    v += [r]
v = v[1:]
v = [list(map(float, line.split(','))) for line in v]


print(len(v))
w = list(map(lambda i : [s[i] for s in v], range(0, len(v))))
ns = w[0]
w = w[1:]
file.close()

ms = [  
        ("1", 0, 1, "#0533ed"),
        ("2",1, 4, "#00FF00"),
        ("3A", 2, 2, "#FF0000"),
        ("3B", 3, 4, "#0000FF"),
        ("4A", 4, 3, "#777777"),
        ("4B", 5, 4, "#FF00FF"),
       #  ("5",  6, 4)  "#00FFFF")
    ] 




plot_eps =  f"""
set terminal png enhanced font 'Verdana,9';
set output 'plot_eps_max_all.png';

set size ratio 0.5;

set format xy "%g";

set logscale xy 10

set grid;
"""
# f(x) = x**(-{p})

# plot """

ps = set([p for (_, _, p, _) in ms])
print(ps)


for p in ps:
    plot_eps += f"f{p}(x) = x**(-{p})\n"

plot_eps += "plot "


for (name, i, p, color) in  ms:


    file = open(f"max_exp_{name}.txt", "w+")
    #file_plt = open(f"max_exp_{name}.plt", "w+") 



    for (n, v) in zip(ns, w[i]):
        file.write(str(n) + "\t" + str(v) + "\n")
   

    plot_eps += f''' "max_exp_{name}.txt" title "eps {name}" with lines linestyle 1 linecolor rgb '{color}', '''
    #file_plt.write(plot_eps + s)

  
    
    file.close()
    #file_plt.close()
for p in p:
    

file_plt = open(f"max_exp_all.plt", "w+") 
file_plt.write(plot_eps)
file_plt.close()
os.system(f"gnuplot max_exp_all.plt")
