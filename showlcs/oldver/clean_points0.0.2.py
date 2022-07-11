#!/usr/bin/python3
#
#   Oleg Ryzhov             oleryz@st.amu.edu.pl
#
#   Clean_Points v. pre-alpha 0.0.2
#   29.05.2022
#

from glob import *
from string import *
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

list_points_glob = []


def show_xy(event):
    # click x-value
    xdata_click = event.xdata
    ydata_click = event.ydata
    global list_points_glob
    list_points_glob.append([xdata_click, ydata_click])
    print(xdata_click, ydata_click)

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

    def __init__(self, inputfile_path, atlfile_in_path, atlfile_out_path, mode, zoom):

        self.mode = mode
        self.zoom = zoom
        self.inputfile_path = inputfile_path
        self.atlfile_in_path = atlfile_in_path
        self.atlfile_out_path = atlfile_out_path
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
        self.NoOfLightcurves = len(self.LightCurves)
        self.Switch = ''

        self.retrunc_part = []  # massive in the future
        self.junk_chunks = []
        self.junk_chunk_len = []
        self.lc_len = []
        self.total_chunks = []
        self.content = []
        self.line_list = []

    def prepare(self):
        if len(sys.argv) == 1:
            self.Switch = ""
        if len(sys.argv) == 2:
            self.Switch = sys.argv[1]

    def reading_atl_file_full(self):

        with open(self.atlfile_in_path) as f:
            self.content = f.readlines()

        self.atl_junk = [item for item in self.content if item[0] != ' ']
        #print(self.atl_junk)
        content = [x.strip() for x in self.content]
        for line in content:
            if line.startswith("OBJECT........:"):
                self.AsteroidName = line.split()[1]
            if line.startswith("LT CORRECTED..:"):
                corr = line.split()[2].strip()
            if line.startswith("OBSERVING TIME:"):
                self.CalDate.append(line.split()[3])
            if line.startswith("ASPECT DATA...:"):
                d_sun = float(line.split()[2])
                d_earth = float(line.split()[3])
                self.PhaseAngle.append(float(line.split()[4]))
                self.R.append(line.split()[2])
                self.Delta.append(line.split()[3])
                self.PA.append(line.split()[4])
                self.Lambda.append(line.split()[5])
                self.Beta.append(line.split()[6])
            if line.startswith("FILTER........:"):
                self.Filter.append(line.split()[1])

            if line.startswith("DATA:"):
                lc = []
                MaxTime = 0
                MinTime = 3000000
                MaxMag = -1000
                MinMag = 1000
            if line.startswith("==============="):
                lc = np.array(lc)
                self.LightCurves.append(lc)
                self.lc_len.append(len(lc))
                self.retrunc_part.append(HelpingModule.retrunc(lc[0][0], -2))
                if self.Ampl < MaxMag-MinMag:
                    self.Ampl = MaxMag-MinMag
                if self.TimeSpan < MaxTime-MinTime:
                    self.TimeSpan = MaxTime-MinTime

            
            tokens = line.split()
            if tokens[0][0] == '!':
                lc.append([0, 0]) # two mods: for plot and for chi-2
                pass
            else:
                JD = float(tokens[0])
                mag = float(tokens[1])
                lc.append([JD, mag])

                if MaxTime < JD:
                    MaxTime = JD
                if MinTime > JD:
                    MinTime = JD

                if MaxMag < mag:
                    MaxMag = mag
                if MinMag > mag:
                    MinMag = mag
        return self.LightCurves
    
    
    def chunk_maker(self):
        
        data_index = [0]
        inds = [ind for ind, item in enumerate(self.atl_junk) if item[:5] == 'DATA:']
        for item in inds:
            data_index.append(item+1)
        i = 1
        #print(data_index)
        for i in range(len(data_index)):
            self.junk_chunks.append(self.atl_junk[data_index[i-1]:data_index[i]])
            self.junk_chunk_len.append(len(self.atl_junk[data_index[i-1]:data_index[i]])) #+1
        self.junk_chunks.pop(0)
        self.junk_chunk_len.pop(0)
        self.junk_chunk_len[0] -= 1

        
        for i in range(len(self.junk_chunk_len)):
            self.total_chunks.append(self.junk_chunk_len[i] + self.lc_len[i])
        print(self.total_chunks)

        #print(self.junk_chunks)
            

    def taking_squares(self):
        global list_points_glob
        list_points = list_points_glob

        #print(list_points)
        if int(self.zoom) == 0:
            pass
        elif int(self.zoom) == 1:
            clean_zoom = [item for ind, item in enumerate(list_points
                                                          ) if ind % (int(self.mode) + int(self.zoom)) != 0]
        else:
            print(
                'Please, insert valid zoom-number (0 or 1 are available)\nAnd follow the updates!')
            quit()

        rect_list = []
        i = 0
        while i < len(clean_zoom):
            rect_list.append(clean_zoom[i:i+int(self.mode)])
            i += int(self.mode)

        for item in rect_list:
            xcoord = []
            ycoord = []
            for item in item:
                xcoord.append(float(item[0]) + self.retrunc_part[0]) #refrence between graphs, points and retrunced part
                ycoord.append(float(item[1]))

            x_min = min(xcoord)
            x_max = max(xcoord)
            y_min = min(ycoord)
            y_max = max(ycoord)
            xcoord.clear()
            ycoord.clear()
            rect_coord = [x_min, x_max, y_min, y_max]
            #print(rect_coord)  
            self.rectangulars.append(rect_coord)

    def find_points(self):

        for i in range(len(self.LightCurves)):

            dict_out = dict()
            coord = dict()
            JD = []
            for item in self.LightCurves[i]:
                JD.append(item[0])
                coord.update({float(item[0]): float(item[1])})
            
            #print(coord)

            for item in self.rectangulars:
                dict_out.update({key: value for (key, value) in coord.items(
                ) if item[0] < key < item[1] and item[2] < value < item[3]})

            indexes = []
            for key in dict_out:
                try:
                    indexes.append(list(coord).index(key) + 1)
                except:
                    'Mistake!'
            
            self.data_out_lc.update({i : indexes})

            coord.clear()
            dict_out.clear()
            JD.clear()

    def get_line(self):
        print(self.data_out_lc)
        for key in self.data_out_lc:
            for item in self.data_out_lc[key]:
                if key == 0:
                    self.line_list.append(self.junk_chunk_len[0] + item)
                elif key == 1:
                    self.line_list.append(self.total_chunks[0] + self.junk_chunk_len[key] + item)
                else:
                    self.line_list.append(sum(self.total_chunks[:key-1]) + self.junk_chunk_len[key] + item) 
            
    def output(self):
        out = [f'!{item[1:]}' if self.line_list.count(ind) > 0 else item for ind, item in enumerate(self.content)]
        #print(self.line_list)
        #print(out)
        with open(self.atlfile_out_path, 'w') as output:
            for item in out:
                output.writelines(f'{item}')  

        #print(self.retrunc_part)
        #print(self.junk_chunks)
        #print(self.junk_chunk_len)
        #print(self.lc_len)

    def plotting(self):
        NoOfLightcurves = len(self.LightCurves)

        grid_size = (5, 5)
        Date = []
        i5 = 0

        print("    Details of lightcurves in the ",
              self.atlfile_in_path, " file")
        print("===============================================================")
        print("No      Date    Filter    r        Delta      PA   lambda  beta")
        print("===============================================================")
        for lc in self.LightCurves:
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

            if NoOfLightcurves == 1:
                cols = 1
                rows = 1
            if NoOfLightcurves == 2:
                cols = 2
                rows = 1
            if NoOfLightcurves == 3:
                cols = 3
                rows = 1
            if NoOfLightcurves == 4:
                cols = 4
                rows = 1
            if NoOfLightcurves > 4:
                cols = 5
                rows = NoOfLightcurves // cols + 1

            if rows == 4:
                figsize = (16, 8)
            if rows == 3:
                figsize = (16, 6)
            if rows == 2:
                figsize = (16, 4)
            if rows == 1:
                figsize = (16, 3)

            if NoOfLightcurves == 1:
                figsize = (4, 3)
            if NoOfLightcurves == 2:
                figsize = (8, 3)
            if NoOfLightcurves == 3:
                figsize = (12, 3)

            gs = gridspec.GridSpec(rows, cols)
            gs.update(hspace=0.4)

            fig1 = plt.figure(num=1, figsize=figsize)
            ax = []

            for lc in self.LightCurves:

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

            plt.suptitle("Lightcurves for %s" %
                         (self.AsteroidName), fontsize=10)
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.9,
                                wspace=None, hspace=None)

            AsteroidPlots = self.atlfile_in_path.split(".")[0]
            plt.savefig(AsteroidPlots+"_rev.pdf")

            print("Close the graphic window to proceed.")
            print("Its content has been saved in " +
                  AsteroidPlots+"_rev.pdf file")

            fig1.canvas.mpl_connect('button_press_event', show_xy)
            plt.show()


if __name__ == '__main__':

    print("Select the ATL file to work with")
    atl_files = glob('*.atl')

    i = 0
    for item in atl_files:
        
        print(f'{i+1}  {atl_files[i]}')
        i += 1

    Ans = input('You choose file number...: ')
    atlfile_in_path = atl_files[int(Ans) - 1]
    # one-point mode smwr in future
    mode = input('How many orientation points do you want to use (4, 3, 2): ')
    # multiple zoom-modes
    zoom = input(
        'How many zooms do you prefer to use (now available only 1): ')
    atlfile_out_path = atlfile_in_path[:-4] + '_mod.atl'
    obj = AutoClean('./input.txt', atlfile_in_path,
                    atlfile_out_path, mode, zoom)
    obj.prepare()
    obj.reading_atl_file_full()
    obj.chunk_maker()
    obj.plotting()
    obj.taking_squares()
    obj.find_points()
    obj.get_line()
    obj.output()
