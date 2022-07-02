#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 12:13:16 2018

@author: Vincent
"""
#%%
from __future__ import division


import os, sys

path = os.getcwd()
sys.path.append(path + "/Calibration_SW/")
sys.path.append("/Users/Vincent/Documents/FireBallIMO/")
from IPython.display import Image as imdisplay


from Calibration.focustest import *  # PlotFocus2DGuider

FUV = Table.read(
    "/Users/Vincent/NextcloudBackUp/FIREBALL/TestsFTS2018-Flight/data/snape/180821/thrufocus/TotalThroughfocus.csv"
)
F1 = Table.read(
    "/Users/Vincent/NextcloudBackUp/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F1_119/TotalThroughfocus.csv"
)
F2 = Table.read(
    "/Users/Vincent/NextcloudBackUp/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F2_-161/TotalThroughfocus.csv"
)
F3 = Table.read(
    "/Users/Vincent/NextcloudBackUp/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F3_-121/TotalThroughfocus.csv"
)
F4 = Table.read(
    "/Users/Vincent/NextcloudBackUp/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F4_159/TotalThroughfocus.csv"
)


def DefineBestActuator(table):
    mean = np.nanmean(
        np.array(
            np.array(
                table[
                    "Best sigma", "Best EE50", "Best EE80", "Best Maxpix", "Best Varpix"
                ]
            ).tolist()
        ),
        axis=1,
    )
    var = np.nanstd(
        np.array(
            np.array(
                table[
                    "Best sigma", "Best EE50", "Best EE80", "Best Maxpix", "Best Varpix"
                ]
            ).tolist()
        ),
        axis=1,
    )
    table["MeanBestActuator"] = mean
    table["VarBestActuator"] = var
    return table


FUV = DefineBestActuator(FUV)
F1 = DefineBestActuator(F1)
F2 = DefineBestActuator(F2)
F3 = DefineBestActuator(F3)
F4 = DefineBestActuator(F4)
#
# table=FUV
# plt.plot()
##plt.plot(table['x'],table['Best EE80'],'x',label = 'EE80')
##plt.plot(table['x'],table['Best EE50'],'+',label = 'EE50')
##plt.plot(table['x'],table['Best sigma'],'.',label = 'sigma')
# plt.errorbar(table['y'],table['MeanBestActuator'],  fmt='o',yerr = table['VarBestActuator'],label = 'Mean')
##plt.plot(table['x'],table['Best Maxpix'],'o',label = 'sigma')
##plt.plot(table['x'],table['Best Varpix'],'o',label = 'sigma')
# plt.legend()
# plt.show()
#
# table = FUV
# X161,Y161,Z161, ax161, Cguider161 = fit_quadratic_curve(table['x'],table['y'],table['MeanBestActuator'],sigma_z = table['VarBestActuator'], n=100,order=1)
# plt.show()


focus = Table.read(
    "/Users/Vincent/Documents/FireBallPipe/Calibration/Focus_180901.csv", format="ascii"
)


def ActuatorAC2Sky(table, column, new_column):
    new_table = table.copy()
    new_table[new_column] = table[column] / 2
    return new_table


def printDiff(tab, sky, ac):
    table = tab.copy()
    table["diff"] = (table[sky] - table[ac]) * 2
    for mask in np.array(range(4)) + 1:
        print("Mask = %i" % (mask))
        print(table[table["mask"] == mask]["star_id", "diff"])
        mean = table[table["mask"] == mask]["diff"].mean()
        std = table[table["mask"] == mask]["diff"].var()
        print("Mean = %0.3f, Std = %0.3f" % (mean, std))
        print("\n")
    return table


def printSum(tab, sky, ac):
    table = tab.copy()
    table["Sum"] = table[sky] + table[ac]
    for mask in np.array(range(4)) + 1:
        print("Mask = %i" % (mask))
        print(table[table["mask"] == mask]["star_id", "Sum"])
        mean = table[table["mask"] == mask]["Sum"].mean()
        std = table[table["mask"] == mask]["Sum"].var()
        print("Mean = %0.3f, Std = %0.3f" % (mean, std))
        print("\n")
    return table


#
# New_focus = ActuatorAC2Sky(focus, 'bestA_Autocol_180825', 'bestA_Autocol_180825/2')
# table = printDiff(New_focus,sky='focus_oc',ac='bestA_Autocol_180825/2')
#
#
# New_focus2 = ActuatorAC2Sky(focus, 'bestA_Autocol_180831', 'bestA_Autocol_180831/2')
# table = printDiff(New_focus2,sky='focus_oc',ac='bestA_Autocol_180831/2')
#
#
# table = printDiff(New_focus,sky='bestA_Autocol_180831',ac='UV_focus_180831')
# table = printDiff(New_focus,sky='UV_focus_180831', ac='bestA_Autocol_180831')
#
# a = printSum(table,sky='diff', ac='focus_oc')
# table = printDiff(New_focus,sky='bestA_Autocol_180825',ac='UV_focus_180825')
#
#

table["UV-Vis_180825_sky"] = (
    table["UV_focus_180825"] - table["bestA_Autocol_180825"]
) * 2
table["UV-Vis_180831_sky"] = (
    table["UV_focus_180831"] - table["bestA_Autocol_180831"]
) * 2

table["UV-Vis_180825"] - table["UV-Vis_180831"]

table["BestUV_from_OC_180825"] = table["focus_oc"] + table["UV-Vis_180825_sky"]
table["BestUV_from_OC_180831"] = table["focus_oc"] + table["UV-Vis_180831_sky"]

table["BestUV_from_OC_180825"] - table["BestUV_from_OC_180831"]


# %%
