#!/usr/bin/python3
#
#   Oleg Ryzhov             oleryz@st.amu.edu.pl
#
#   Show LightCurves v. pre-alpha 1.8
#   29.05.2022
#

try:
    from glob import *
except ImportError:
    print('Import error. Please, write following in command line: pip install glob')

try:
    from string import *
except ImportError:
    print('Import error. Please, write following in command line: pip install string')

try:
    import numpy as np
except ImportError:
    print('Import error. Please, write following in command line: pip install numpy')

try:
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
except ImportError:
    print('Import error. Please, write following in command line: pip install matplotlib')

try:
    import re
except ImportError:
    print('Import error. Please, write following in command line: pip install re')

try:
    import random
except ImportError:
    print('Import error. Please, write following in command line: pip install random')

try:
    import os
except ImportError:
    print('Import error. Please, write following in command line: pip install os')

import sys

list_points_glob = []
t = 0


def show_xy(event):
    global list_points_glob
    if event.button == 3:
        print('Clicked the RB')
        xdata_click = event.xdata
        ydata_click = event.ydata
        list_points_glob.append([xdata_click, ydata_click])
        print(xdata_click, ydata_click)
    if event.button == 1:
        pass


class HelpingModule:

    def trunc(JD, decimals):
        multiplier = 10**decimals
        rounded_JD = int(JD * multiplier) / multiplier
        MJD = JD - rounded_JD
        return round(MJD, 6)

    def retrunc(JD, decimals):
        multiplier = 10**decimals
        rounded_JD = int(JD * multiplier) / multiplier
        return rounded_JD


class AutoClean(HelpingModule):

    def __init__(self, DataFile, cache_dir):

        self.DataFile = DataFile
        self.atl_outpath = DataFile[:-4] + '@.atl'
        self.pdf_file = []
        self.data_out_lc = {}
        self.rectangulars = []
        self.atl_junk = []
        self.AsteroidName = ''
        self.docheck = 0
        self.Filter = []
        self.PhaseAngle = []
        self.CalDate = []
        self.AspectData = []
        self.R = []
        self.Delta = []
        self.PA = []
        self.Lambda = []
        self.Beta = []
        self.Ampl = 0.0
        self.TimeSpan = 0.0
        self.LightCurves = []
        self.NoOfLightcurves = 0
        self.Switch = ''
        self.retrunc_part = []
        self.junk_chunks = []
        self.junk_chunk_len = []
        self.lc_len = []
        self.total_chunks = []
        self.content = []
        self.line_list = []
        self.exception = []
        
        self.zerotime = []
        self.content = []
        self.content_atl = []
        self.cap = []
        self.lc_dict = {}
        self.mags_list = []
        self.phases_list = []
        self.ind_list = []

        self.ind_coord = {}
        self.ind_out = {}
        self.index_dict_out = {}

        self.lc_str = []
        self.lc_str_em = []
        self.lc_str_all = []
        self.cache_dir = cache_dir

    def prepare(self):
        if len(sys.argv) == 1:
            self.Switch = ""
        if len(sys.argv) == 2:
            self.Switch = sys.argv[1]

    def read_atl_lightcurves(self):
        
        with open(self.DataFile) as f:
            self.content = [item.strip() for item in f.readlines()]
            lcnjunk = []
            for item in self.content:
                if item.startswith('=') == False:
                    lcnjunk.append(item)
                else:
                    self.content_atl.append(lcnjunk)
                    lcnjunk = []
            i = 1
            for item in self.content_atl:
                t = -1
                for line in item[13:25]:
                    if line.startswith('INFORMATION...:'):
                        t += 1
                self.cap.append(item[:19+t])
                mid_lc = [' ' + item if item.startswith('!') == False else item for item in item[19+t:]]
                self.lc_str.append(mid_lc)
                mid_lc = [ind for ind, item in enumerate(item[19+t:]) if item.startswith('!') == True]
                self.lc_str_em.append(mid_lc)
                mid_lc = [' ' + item if item.startswith('!') == False else item for item in item[19+t:]]
                self.lc_str_all.append(mid_lc)

                self.AsteroidName = item[0].split()[1]
                self.Filter.append(item[8].split()[1])
                self.RelPlot = item[10].split()[2].strip()
                self.CalDate.append(item[16+t].split()[3])
                self.R.append(item[17+t].split()[2])
                self.Delta.append(item[17+t].split()[3])
                self.PA.append(float(item[17+t].split()[4]))
                self.Lambda.append(item[17+t].split()[5])
                self.Beta.append(item[17+t].split()[6])
                lc = []
                MaxTime = 0
                MinTime = 3000000
                MaxMag = -1000
                MinMag = 1000
                for item in item[19+t:]:
                    tokens = item.split()
                    if tokens[0][0] == '!':
                        pass
                    else:
                        try: 
                            JD = float(tokens[0])
                            mag = float(tokens[1])
                            magsig = float(tokens[2])
                        except:
                            print(f'Escaping the program. There is a mistake in {self.DataFile}, correct it: {item}')
                            exit()
                        lc.append([JD, mag, magsig, i])
                        i += 1
                        if MaxTime < JD:
                            MaxTime = JD
                        if MinTime > JD:
                            MinTime = JD

                        if MaxMag < mag:
                            MaxMag = mag
                        if MinMag > mag:
                            MinMag = mag
                self.LightCurves.append(np.array(lc))
                self.lc_len.append(len(lc))
                self.junk_chunk_len.append(19+t)
                self.total_chunks.append(len(item))
                if self.Ampl < MaxMag-MinMag:
                    self.Ampl = MaxMag-MinMag
                if self.TimeSpan < MaxTime-MinTime:
                    self.TimeSpan = MaxTime-MinTime
            for item in self.LightCurves:
                self.zerotime.append(item[0,0])

    def taking_squares(self):
        global list_points_glob
        list_points = list_points_glob
        rect_list = []
        i = 0
        
        while i < len(list_points):
            rect_list.append(list_points[i:i+2])
            i += 2
            
        for item in rect_list:
            xcoord = []
            ycoord = []
            for item in item:
                xcoord.append(float(item[0])) #refrence between graphs, points and retrunced part
                ycoord.append(float(item[1]))

            x_min = min(xcoord)
            x_max = max(xcoord)
            y_min = min(ycoord)
            y_max = max(ycoord)
            xcoord.clear()
            ycoord.clear()  
            self.rectangulars.append([x_min, x_max, y_min, y_max])

    def transfer_data(self):
        for i in range(len(self.ind_list)):
            phases = []
            mags2 = []
            ind = []
            phases.append(self.phases_list[i].tolist())
            mags2.append(self.mags_list[i].tolist())
            ind.append(self.ind_list[i].tolist())
            for j in range(len(phases[0])):
                self.ind_coord.update({ind[0][j] : [phases[0][j], mags2[0][j], i]})

    def find_points(self):

        for item in self.rectangulars:
            self.index_dict_out.update({key : value[2] for (key,value) in self.ind_coord.items() if item[0] < value[0] < item[1] and item[2] < value[1] < item[3]})
        
        self.index_dict_out = {key : value for key, value in sorted(self.index_dict_out.items(), key = lambda item: item[1])}
        if len(self.index_dict_out) == 0:
            return
        for i in self.index_dict_out:
            self.ind_out.setdefault(self.index_dict_out[i], []).append(i)

    def get_line(self):
        for key in self.ind_out:
            line_list = []
            for item in self.ind_out[key]:
                if key == 0:
                    line_list.append(item - 1)
                elif key == 1:
                    line_list.append(item - 1 - self.lc_len[0]) 
                else:
                    line_list.append(item - 1 - sum(self.lc_len[:key]))
            self.lc_dict.update({key : line_list})
        print(self.lc_dict)
        
    def changing_str(self):
        em = []
        for i in range(len(self.lc_str_em)):
            em.append([item for ind, item in enumerate(self.lc_str_em[i])])
        for key in self.lc_dict:
            for item in self.lc_dict[key]:
                ind = int(item)
                for i in range(len(em[key])):
                    if ind >= em[key][i]:
                        ind += 1
                    else:
                        break
                try:
                    self.lc_str_all[key][ind] = '!' + self.lc_str_all[key][ind][1:]
                except:
                    print('Some indexing problem')

    def rewriting_atl(self):
        with open(self.atl_outpath, 'w') as out_atl:
            for i in range(len(self.cap)):
                for item in self.cap[i]:
                    out_atl.writelines(f'{item}\n')
                for item in self.lc_str_all[i]:
                    out_atl.writelines(f'{item}\n')
                out_atl.writelines('================\n')
            out_atl.writelines(' END OF FILE\n')


    def plotting(self):
        self.NoOfLightcurves = len(self.LightCurves)
        splitted_lc = [self.LightCurves[i: i + 19] for i in range(0, self.NoOfLightcurves, 19)]
        for item in splitted_lc:
            len_lightcurves = len(item)
            print(len_lightcurves)
            grid_size = (5, 5)
            Date = []
            i5 = 0

            print("    Details of lightcurves in the ",
                self.DataFile, " file")
            print("===============================================================")
            print("No      Date    Filter    r        Delta      PA   lambda  beta")
            print("===============================================================")
            for lc in item:
                # Date.append(str(jd2gcal(0.0, float(lc[0,0])[0:3]))
                Y, M, D = self.CalDate[i5].split('-')
                Y1 = int(Y)
                M1 = int(M)
                D1 = int(D)
                Y1 = int(Y)
                M1 = int(M)
                D1 = int(D)
                PA1 = float(self.PA[i5])
                Lambda1 = float(self.Lambda[i5])
                Beta1 = float(self.Beta[i5])
                print("%2i   %s   %s      %s     %s   %4.1f  %5.1f  %4.1f" %
                    (i5+1, self.CalDate[i5], self.Filter[i5], self.R[i5], self.Delta[i5], PA1, Lambda1, Beta1))

                i5 = i5+1
            print("===============================================================")

            if self.Switch != "-n":
                # -------------------------------------------------
                # 4.  Plot lightcurves
                # -------------------------------------------------

                # define the figure size and grid layout properties
                # keep the number of columns fixed at 5, adjust figsize
                #
                i5 = 0

                if len_lightcurves == 1:
                    cols = 1
                    rows = 1
                if len_lightcurves == 2:
                    cols = 2
                    rows = 1
                if len_lightcurves == 3:
                    cols = 3
                    rows = 1
                if len_lightcurves == 4:
                    cols = 4
                    rows = 1
                if len_lightcurves > 4:
                    cols = 5
                    rows = len_lightcurves // cols + 1
                if rows == 4:
                    figsize = (16, 8)
                if rows == 3:
                    figsize = (16, 6)
                if rows == 2:
                    figsize = (16, 4)
                if rows == 1:
                    figsize = (16, 3)

                if len_lightcurves == 1:
                    figsize = (4, 3)
                if len_lightcurves == 2:
                    figsize = (8, 3)
                if len_lightcurves == 3:
                    figsize = (12, 3)

                gs = gridspec.GridSpec(rows, cols)
                gs.update(hspace=0.4)

                fig1 = plt.figure(num=1, figsize=figsize)
                ax = []

                for lc in item:

                    col = np.remainder(i5, 5)
                    row = np.floor_divide(i5, 5)
                #    if NoOfLightcurves==1:
                #      col=1
                #      row=1

                    XMin = lc[0, 0]-np.trunc(lc[0, 0]/100.0)*100.0
                    XMax = XMin+self.TimeSpan
                    AveMag = (np.amax(lc[:, 1])+np.amin(lc[:, 1]))/2.0
                    YMax = AveMag+0.6*self.Ampl
                    YMin = AveMag-0.6*self.Ampl
                    ax.append(fig1.add_subplot(gs[row, col]))
                    # ax[-1].get_xaxis().set_visible(False)
                    ax[-1].set_xlim([XMin, XMax])
                    ax[-1].set_ylim([YMin, YMax])
                    # ax[-1].get_yaxis().set_fontsize(20)
                    plt.gca().invert_yaxis()
                    plt.tick_params(labelsize=8)
                    ax[-1].set_title("%3i    %s %s %3.1f%s %s" % (i5+1, self.CalDate[i5],
                                    "PA=", float(self.PA[i5]), ",", self.Filter[i5]), fontsize=8)
                    ax[-1].plot(lc[:, 0]-np.trunc(lc[:, 0]/100.0)
                                * 100.0, lc[:, 1], 'ro')
                    i5 = i5+1
                    phase=lc[:, 0]-np.trunc(lc[:, 0]/100.0)*100.0
                    mags= lc[:,1]
                    ind=lc[:,3]

                    self.mags_list.append(mags)
                    self.phases_list.append(phase)
                    self.ind_list.append(ind)

                plt.suptitle("Lightcurves for %s" %
                            (self.AsteroidName), fontsize=10)
                plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.9,
                                    wspace=None, hspace=None)

                AsteroidPlots = self.DataFile.split(".")[0]
                pdf = AsteroidPlots + randomising() + '_rev.pdf'
                self.pdf_file.append(pdf)
                plt.savefig(pdf)

                print("Close the graphic window to proceed.")
                print(f'Its content has been saved in {pdf} file')

                fig1.canvas.mpl_connect('button_press_event', show_xy)
                plt.show()
        
    def caching(self):
        cur_dir = os.getcwd()
        atlfile_inpath = os.path.join(cur_dir, self.DataFile)
        atlfile_outpath = os.path.join(self.cache_dir, self.DataFile)
        os.rename(atlfile_inpath, atlfile_outpath)

        for pdf in self.pdf_file:
            pdffile_inpath = os.path.join(cur_dir, pdf)
            pdffile_outpath = os.path.join(cache_dir, pdf)
            os.rename(pdffile_inpath, pdffile_outpath)
        print('Caching is finished :3')

def randomising():
    random.seed()
    numbers = np.arange(0, 9, 1)
    name = ''
    for i in range(0, 10):
        name = name + str(random.choice(numbers))
    return name

def clean_cache(cache_dir, num_sym, cache_answer): 
    files_out = [item for item in os.listdir(cache_dir) if len(item) == num_sym or len(re.sub(r'_rev', '', item)) + 10 == num_sym]
    files_del = [item for item in os.listdir(cache_dir) if len(item) != num_sym and len(re.sub(r'_rev', '', item)) + 10 != num_sym]
    for file in files_out:
        os.rename(os.path.join(cache_dir, file), os.path.join(os.getcwd(), file))
    if cache_answer == 'N':
        user_dir = str(os.path.join(os.getcwd(), randomising()))
        os.mkdir(user_dir)
        for file in files_del:
            try:
                os.rename(os.path.join(cache_dir, file), os.path.join(user_dir, file))
            except:
                print('Passing...')
                pass
    for file in os.listdir(cache_dir):
        os.remove(os.path.join(cache_dir, file)) 
    out_names = randomising() 
    for file in os.listdir(os.getcwd()):
        signed = re.search('@', file)
        if signed:
            print(file)
            file_mod = out_names + file[-4:]
            os.rename(os.path.join(os.getcwd(), file), os.path.join(os.getcwd(), file_mod))
            if file_mod[-4:] == '.atl':
                show_blockchain(file_mod, cache_dir)

def blockchain(DataFile, cache_dir):
    obj = AutoClean(DataFile, cache_dir)
    obj.prepare()
    obj.read_atl_lightcurves()
    obj.plotting()
    obj.transfer_data()
    obj.taking_squares()
    obj.find_points()
    obj.get_line()
    obj.changing_str()
    obj.rewriting_atl()
    obj.caching()

def show_blockchain(DataFile, cache_dir):
    obj = AutoClean(DataFile, cache_dir)
    obj.prepare()
    obj.read_atl_lightcurves()
    obj.plotting()

def recall():
    answer = input('Do you want to recall the program. [Y/N] ')
    if answer == 'Y':
        return True
    return False


if __name__ == '__main__':

    print("Select the ATL file to work with")
    atl_files = glob('*.atl')

    i = 0
    for item in atl_files:
        print(f'{i+1}  {atl_files[i]}')
        i += 1

    Ans = input('You choose file number...: ')
    DataFile = atl_files[int(Ans) - 1]
    num_sym = len(DataFile)
    cur_dir = os.getcwd()
    cache_dir = os.path.join(cur_dir, DataFile[:-4] + '_cache')
    try:
        os.mkdir(cache_dir)
    except:
        pass
    blockchain(DataFile, cache_dir)
    answer = recall()
    while answer != False:
        DataFile = DataFile[:-4] + '@.atl'
        blockchain(DataFile, cache_dir)
        answer = recall()
    
    cache_answer = input('Do you want cache-file to be cleared? [Y/N] ')

    clean_cache(cache_dir, num_sym, cache_answer)
    print('Thanks for using the showlcs program <3')
    
