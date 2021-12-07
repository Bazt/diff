import os

ns = [10, 80, 1280]
ms = [  
        ("1", "#0533ed"),
        ("2", "#00FF00"),
        #("3A", "#FF0000"),
       # ("3B", "#0000FF"),
        ("4A", "#FFFF00"),
        ("4B", "#FF00FF"),
       # ("5",  "#00FFFF")
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

