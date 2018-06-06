## EZ 2018 Only print statement updates were needed to make this file Python 3 compatible. earl.io
'''
Copyright (C) 2011, 2014 by Xuebing Wu and Robert C. Berwick

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

GUI
'''

import os,sys,wx,matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

model = sys.argv[1]

if model == 'coalesce':
    from coalesce_with_mutation import *
    title = 'Coalesce with mutation'
    paraName = ['population size (N)','num of lineages','mutation rate']
    paraValue = [25,6,0.01]
    paraType = ['int','int','float']
        
elif model == 'sex':
    from sexsel import *
    title = 'Sexual selection'
    paraName = ['population size (N)','num of generations','mate selection preference','survival probability','male genotype number','female genotype number','fix population size',]
    paraValue = [1000,15,'0.25\t0.25\t0.25\t0.25\n0\t0\t0.5\t0.5\n0.25\t0.25\t0.25\t0.25\n0\t0\t0.5\t0.5',1.0,'250,250,250,250','250,250,250,250',True]
    paraType = ['int','int','multipleline','float','string','string','bool']    
        
elif model == 'drift':
    from diploid_simulation import *
    title = 'drift, mutation, and migration'
    paraName = ['population size (N)',
            'num of populations',
            'num of generations',
            'starting allele A freq',
            'mutation rate (A2a,a2A)',
            'relative fitness (AA,Aa,aa)',
            'frac. migrants each generation',
            'allele A freq in migrants',
            'frequency-dependent selection (s)']
    paraValue = [10,20,100,0.5,'0,0','1,1,1','0.0','1.0','-1']
    paraType = ['int','int','int','float','string','string','float','float','float']    

elif model == 'muller':
    from muller1 import *
    title = 'muller ratchet'
    paraName = ['population size (N)', 'U', 'fitness','mutmax']
    paraValue = [100,0.1,0.9,30]
    paraType = ['int','float','float','int']
    
def run(fig,para,model):
    if model == 'coalesce':        
        lineage,coalesce_time,mutation_time,coalesce_event,mutation_event = coalesce_with_mutation(para[0],para[1],para[2])
        drawtree(fig,coalesce_time,mutation_time,coalesce_event,mutation_event)
    elif model =='sex':
        male_genotype_freqs,female_genotype_freqs, population_size = simulation(para[0],para[1],para[2],para[3],para[4],para[5],para[6]=='yes')
        draw(fig,male_genotype_freqs,female_genotype_freqs)
    elif model == 'drift':
        freq_matrix = simulation(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8])
        draw(fig,freq_matrix)
    elif model == 'muller':
        generations,minimum = muller(para[0],para[1],para[2],para[3])
        axes = fig.add_subplot(111)
        axes.plot(generations,minimum,'+')
        axes.set_xlim(0,max(generations))
        axes.set_xlabel('generation')
        axes.set_ylabel('minimum # of mutations')          

            
# below are code for GUI                
class BarsFrame(wx.Frame):
    ''' 
    The main frame of the application
    '''
    title = title
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)
        
        self.para = paraValue      
        
        self.create_menu()
        self.create_main_panel()
        
        for i in range(len(paraName)):
            #print paraName[i],paraType[i],self.para[i]
            if paraType[i] == 'bool':
                self.textbox[i].SetValue(self.para[i])
            else:
                self.textbox[i].SetValue(str(self.para[i]))            
        self.run_and_draw()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        ''' 
        Creates the main panel with all the controls on it
        '''
        self.panel = wx.Panel(self) 
        
        self.dpi = 100
        self.fig = Figure((7.0, 6.0),facecolor="white", dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        

        self.label = [None]*len(paraName)
        self.textbox = [None]*len(paraName)        
        for i in range(len(paraName)):
            self.label[i] = wx.StaticText(self.panel, -1, paraName[i])
            if paraType[i] == "multipleline":
                self.textbox[i] = wx.TextCtrl(self.panel,size=(200,70),style=wx.TE_MULTILINE)
            elif paraType[i] == 'bool':
                self.textbox[i] = wx.CheckBox(self.panel,size=(200,-1))
            else: 
                self.textbox[i] = wx.TextCtrl(self.panel,size=(200,-1))
            #self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox[i])  
            
        self.drawbutton = wx.Button(self.panel, -1, "RUN")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.hbox.AddSpacer(10)
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        for i in range(len(paraName)):
            self.vbox.Add(self.label[i],0, border=3, flag=flags)
            self.vbox.Add(self.textbox[i], 0, border=3, flag=flags)
        self.vbox.Add(self.drawbutton, 0, border=3, flag=flags)
        self.vbox.AddSpacer(30)
        
        self.hbox.Add(self.vbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.panel.SetSizer(self.hbox)
        self.hbox.Fit(self)

    def run_and_draw(self):
        ''' 
        update
        '''
        print('---------new simulation---------')
        for i in range(len(paraName)):
            if paraType[i] == 'int': 
                self.para[i] = int(self.textbox[i].GetValue())
            elif paraType[i] == 'float':           
                self.para[i] = float(self.textbox[i].GetValue())
            else:
                self.para[i] = self.textbox[i].GetValue()
            print(paraName[i],'(',paraType[i],')',self.para[i])

        self.fig.clear()
        
        run(self.fig,self.para,model)
        
        self.canvas.draw()
    
    def on_draw_button(self, event):
        self.run_and_draw()  
    
    def on_text_enter(self, event):
        self.run_and_draw()

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        
        dlg = wx.FileDialog(
            self, 
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
        
    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = '''
        Population genetics simulation
        Author: Xuebing Wu (wuxbl@mit.edu)
        '''
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
       
        
if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = BarsFrame()
    app.frame.Show()
    app.MainLoop()

