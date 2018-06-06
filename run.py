#!/usr/bin/python
## EZ 2018 this file required no modifications to work with Python 3. earl.io
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
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
import os,wx

APP_SIZE_X = 300
APP_SIZE_Y = 200

class Menu(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(APP_SIZE_X, APP_SIZE_Y))

        self.panel = wx.Panel(self)
        
        self.drift = wx.Button(self.panel, -1, "Drift, mutation, migration, selection")
        self.Bind(wx.EVT_BUTTON, self.on_drift, self.drift)
        
        self.coal = wx.Button(self.panel, -1, "Coalescence with mutation")
        self.Bind(wx.EVT_BUTTON, self.on_coal, self.coal)

        self.sex = wx.Button(self.panel, -1, "Sexual selection")
        self.Bind(wx.EVT_BUTTON, self.on_sex, self.sex)

        self.muller = wx.Button(self.panel, -1, "Muller's ratchet")
        self.Bind(wx.EVT_BUTTON, self.on_muller, self.muller)
                        

        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.drift,1, wx.LEFT | wx.TOP | wx.GROW)
        box.Add(self.coal,1, wx.LEFT | wx.TOP | wx.GROW )
        box.Add(self.sex, 1, wx.LEFT | wx.TOP | wx.GROW )
        box.Add(self.muller, 1, wx.LEFT | wx.TOP | wx.GROW )                
        self.panel.SetSizer(box)
        
        self.Centre()
        self.ShowModal()
        self.Destroy()

    def on_drift(self,event):
        os.system('python GUI.py drift')
    def on_coal(self,event):
        os.system('python GUI.py coalesce')
    def on_sex(self,event):
        os.system('python GUI.py sex')
    def on_muller(self,event):
        os.system('python GUI.py muller')                        
        
    def OnClose(self, event):
        self.Close(True)


app = wx.App(0)
Menu(None, -1, '7.33 Simulation')
app.MainLoop()
