# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 11:25:20 2018

@author: thomas.s.gresock
"""
import math
import pandas as pd
import numpy as np
import calendar
import time
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import Button
from matplotlib.widgets import MultiCursor
from matplotlib.widgets import RadioButtons
from matplotlib.widgets import TextBox
import matplotlib.dates as mdates
from tkinter import Tk
from tkinter.filedialog import askopenfilenames



csv_name = 'EPOCH_' + str(calendar.timegm(time.gmtime())) + '.csv'
meta_frame = {}
label = 'Files_One'
hz_flag = 'Ten'
rate = 10
hz_dict = {'Ten':{'value':0,'color':'blue','percent':10},
           'Forty':{'value':0,'color':'green','percent':40},
           'Eighty':{'value':0,'color':'red','percent':80},
           'One_hundred':{'value':0,'color':'black','percent':100}}


fig = plt.figure(figsize=(9,7))
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.3)
ax = fig.add_subplot(312)
ax2 = fig.add_subplot(311)
ax3 = fig.add_subplot(313)
ax.set_xlabel('Distance')
ax.set_ylabel('Comp')
ax2.set_xlabel('Distance')
ax2.set_ylabel('dBm')
ax3.set_xlabel('Time')
ax3.set_ylabel('Range')


def vincenty_inverse(point1, point2, miles=False):
    a = 6378137  # meters
    f = 1 / 298.257223563
    b = 6356752.314245  # meters; b = (1 - f)a

    MILES_PER_KILOMETER = 0.621371

    MAX_ITERATIONS = 200
    CONVERGENCE_THRESHOLD = 1e-12
    """
    Extracted from Vincenty module at 
    """

    # short-circuit coincident points
    if point1[0] == point2[0] and point1[1] == point2[1]:
        return 0.0

    U1 = math.atan((1 - f) * math.tan(math.radians(point1[0])))
    U2 = math.atan((1 - f) * math.tan(math.radians(point2[0])))
    L = math.radians(point2[1] - point1[1])
    Lambda = L

    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)

    for iteration in range(MAX_ITERATIONS):
        sinLambda = math.sin(Lambda)
        cosLambda = math.cos(Lambda)
        sinSigma = math.sqrt((cosU2 * sinLambda) ** 2 +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
        if sinSigma == 0:
            return 0.0  # coincident points
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        sigma = math.atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha ** 2
        try:
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        except ZeroDivisionError:
            cos2SigmaM = 0
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LambdaPrev = Lambda
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                               (cos2SigmaM + C * cosSigma *
                                                (-1 + 2 * cos2SigmaM ** 2)))
        if abs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
            break  # successful convergence
    else:
        return None  # failure to converge

    uSq = cosSqAlpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                 (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                 (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
    s = b * A * (sigma - deltaSigma)

    s /= 1000  # meters to kilometers
    if miles:
        s *= MILES_PER_KILOMETER  # kilometers to miles

    return round(s, 6)


def set_hz_flag(label):
    global hz_flag
    if label == '10%':
        hz_flag = 'Ten'
    elif label == '40%':
        hz_flag = 'Forty'
    elif label == '80%':
        hz_flag = 'Eighty'
    elif label == '100%':
        hz_flag = 'One_hundred'
#radio.on_clicked(set_hz_flag)


def choose_transmitter(event):
    ax.clear()
    ax2.clear()
    ax3.clear()
    Tk().withdraw()
    global transmitter_filename
    global receiver_filename
    global transmitter_gps_filename
    global receiver_gps_filename
    global filenames_list
    global trans_mac
    global df
    global left 
    global right
    filenames_list = askopenfilenames()
    for x in filenames_list:
        if x.split('.')[-1] == 'cw14tx':
            transmitter_filename = x
        elif x.split('.')[-1] == 'cw14rx':
            receiver_filename = x
        elif x.split('.')[-1] == 'gpstx':
            transmitter_gps_filename = x
        elif x.split('.')[-1] == 'gpsrx':
            receiver_gps_filename = x
        else:
            print(x)
    tr = transmitter_filename.split('_')[-1].split('.')[0]
    if tr == 'GVW':
        trans_mac = '04:e5:48:01:9b:e0'
        print(tr)
    elif tr == 'WVW':
        trans_mac = '04:e5:48:01:9c:10'
    elif tr == 'RSU':
        trans_mac = '04:e5:48:01:50:d0'
    elif tr == 'MYS':
        trans_mac = '04:e5:48:01:3b:a4'
    df = build_df(transmitter_filename,
                  receiver_filename, 
                  transmitter_gps_filename,
                  receiver_gps_filename,
                  trans_mac)
    left = df.index[0]
    right = df.index[-1]
    ax3.plot(df.index,df.range,marker='.')
    ax3.set_xlim(left,right)
    return

def read_files(transmitter, reciever, gps1, gps2, trans_mac):
    transmit = pd.read_csv(transmitter, header=None,
                           names=[
                                  'Transmit_Time_1',
                                  'ID_Time_1',
                                  'Index_1',
                                  'lat_1',
                                  'lon_1',
                                  'u1_1',
                                  'u2_1',
                                  'msg_length_1',
                                  'u3_1',
                                  'u4_1'                                  
                                  ])
    recieve = pd.read_csv(reciever, header=None,
                          names=[
                                  'Rec_Time_1',
                                  'Sent_Time_1',
                                  'Rec_lat_1',
                                  'Rec_lon_1',
                                  'Rec_heading_1',
                                  'u1_1',
                                  'msg_len_1',
                                  'transmitter_ID_1',
                                  'Transmit_Index_1',
                                  'transmitter_lat_1',
                                  'transmitter_lon_1',
                                  'u2_1',
                                  'u3_1',
                                  'RSSI_1',
                                  'RSSI2_1',
                                  'RSSI3_1',
                                  'RSSI4_1',
                                  'u4_1',
                                  'u5_1'
                                  ],skiprows=2)
    recieve.dropna(inplace=True)
    recieve = recieve[recieve.u5_1 == trans_mac]



    
    gps1 = pd.read_csv(gps1, header=None,
                      names=[
                             'gps_time1',
                             'drop1',
                             'gps_lat1',
                             'gps_lon1',
                             'drop11',
                             'drop21'
                             ])
    gps2 = pd.read_csv(gps2, header=None,
                      names=[
                             'gps_time2',
                             'drop2',
                             'gps_lat2',
                             'gps_lon2',
                             'drop12',
                             'drop22'
                             ])
    return transmit, recieve, gps1, gps2
    
    
def build_df(transmitter, reciever, rtk1, rtk2, trans_mac):
    one, two, three, four = read_files(transmitter, reciever, rtk1, rtk2, trans_mac)

    one = one[
               [
                'Transmit_Time_1',
                'ID_Time_1',
                'Index_1',
                'lat_1',
                'lon_1',
                ]
             ]
               
    two = two[
              [
               'Rec_Time_1',
               'Sent_Time_1',
               'Rec_lat_1',
               'Rec_lon_1',
               'transmitter_ID_1',
               'Transmit_Index_1',
               'transmitter_lat_1',
               'transmitter_lon_1',
               'RSSI_1',
               'RSSI2_1',
              ]
             ]
    two['comp'] = 1

    one = one.set_index('Index_1', drop=False)
    two = two.set_index('Transmit_Index_1', drop=False)

    inter = one.join(two,how='left')

    inter['comp'].fillna(0, inplace=True)
    inter2 = inter
    inter2['Time_in_transit'] = inter2.Rec_Time_1 - inter2.Sent_Time_1
    inter2['Time_in_q'] = inter2.Sent_Time_1 - inter2.Transmit_Time_1
    inter2['Time_from_q_to_rec'] = inter2.Rec_Time_1 - inter2.Transmit_Time_1
    inter = inter.set_index('Transmit_Time_1', drop=False)
    
    three = three.set_index('gps_time1',drop=False)
    four = four.set_index('gps_time2',drop=False)
    
    three = three[['gps_time1','gps_lat1','gps_lon1']]
    four = four[['gps_time2','gps_lat2','gps_lon2']]
    
    
    
    gps_final = three.join(four,how='outer')
    #gps_final = gps_final.join(cb_final,how='outer')

    gps_final = gps_final.join(inter,how='outer')
    gps_final = gps_final[['gps_time1','gps_time2','Transmit_Time_1','Rec_Time_1','Sent_Time_1',
                           'Index_1','Transmit_Index_1','RSSI_1','RSSI2_1','comp',
                           'gps_lat1','gps_lat2','gps_lon1','gps_lon2','Time_in_transit',
                           'Time_in_q','Time_from_q_to_rec']]

    
    gps_final['gps_lat1'] = gps_final.gps_lat1.interpolate()
    gps_final['gps_lat2'] = gps_final.gps_lat2.interpolate()
    gps_final['gps_lon1'] = gps_final.gps_lon1.interpolate()
    gps_final['gps_lon2'] = gps_final.gps_lon2.interpolate()
    
    gps_final['range'] = np.sqrt((6371000.0*(gps_final.gps_lat1 - gps_final.gps_lat2)*0.017453)**2 +
                         (6371000.0*np.cos(gps_final.gps_lat1*0.017435)*((gps_final.gps_lon1 -
                          gps_final.gps_lon2)*0.017453))**2)

    gps_final = gps_final.set_index(pd.to_datetime(gps_final.index, unit='s'))

    return gps_final#, one, inter
    
 









def run1(event):
    ax.clear()
    ax2.clear()
    print(transmitter_filename)
    print(receiver_filename)
    print(transmitter_gps_filename)
    print(receiver_gps_filename)
    '''
    df = build_df(transmitter_filename,
                  receiver_filename, 
                  transmitter_gps_filename,
                  receiver_gps_filename)
    '''
    sub_df = df[left:right]
    df_comp = sub_df[['range','comp']].dropna()
    df_comp_rolling = df_comp.rolling(rate).mean()
    df_RSSI = sub_df[['range','RSSI_1']].dropna()
    ax.plot(df_comp.range,df_comp.comp,marker='.')
    ax.plot(df_comp_rolling.range,df_comp_rolling.comp,color='red',linewidth=1)
    ax2.plot(df_RSSI.range,df_RSSI.RSSI_1,marker='.')
    #ax3.plot(df.index,df.range,marker='.')
    
    
    #print(df)
    return df
    

def onselect(xmin, xmax):
    global left
    global right
    global span_1
    print('start: {} , {}'.format(xmin, xmax))
    left = pd.to_datetime(mdates.num2epoch(xmin),unit='s')
    right = pd.to_datetime(mdates.num2epoch(xmax),unit='s')
    try:
        span_1.remove()
        span_1 = ax3.axvspan(xmin,xmax,facecolor='red',alpha=0.4)
    except:
        span_1 = ax3.axvspan(xmin,xmax,facecolor='red',alpha=0.4)            
    print('start: {} , {}'.format(left, right))
    #print('endpoint: {}, {}'.format(erelease.xdata,erelease.ydata))
    
def submit_rolling(text):
    global rate
    rate = int(text)
    print('rolling mean set to {} datapoints'.format(rate))
    
def onclick_ax(event):
    global hz_dict
    global v_line_40
    global v_line_60
    global v_line_80
    global v_line_100
    #global ix
    if event.inaxes == ax:
        print(event.xdata)
        if hz_flag == 'Ten':
            try:                
                v_line_10.remove()
                v_line_10 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)
            except:
                v_line_10 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)
        elif hz_flag == 'Forty':
            try:                
                v_line_40.remove()
                v_line_40 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)
            except:
                v_line_40 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)
        elif hz_flag == 'Eighty':
            try:
                v_line_80.remove()
                v_line_80 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)
            except:
                v_line_80 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)        
        elif hz_flag == 'One_hundred':
            try:
                v_line_100.remove()
                v_line_100 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)
            except:
                v_line_100 = ax.axvline(event.xdata,color=hz_dict[hz_flag]['color'],linewidth=4)
        else:
            print('broken')
        hz_dict[hz_flag]['value'] = event.xdata
        
        
        
def change_label(text):
    global label
    label = text
    
def move_to_frame_func(event):
    global label
    global meta_frame
    current_sub_frame = {10:hz_dict['Ten']['value'],
                         40:hz_dict['Forty']['value'],
                         80:hz_dict['Eighty']['value'],
                         100:hz_dict['One_hundred']['value']}
    meta_frame[label] = current_sub_frame
    print('{}'.format(meta_frame))
    
    
    
def set_csv_name(csv_text):
    global csv_name
    csv_name = csv_text + '.csv'
    
def write_csv_func(event):
    meta_df = pd.DataFrame(meta_frame)
    meta_df.sort_index()
    meta_df.to_csv(csv_name)
    print('FILE WRITTEN')
    
def close(event):
    plt.close('all')
    quit


cid = fig.canvas.mpl_connect('button_press_event',onclick_ax)

'''RADIO BUTTONS'''
axcolor = 'lightgoldenrodyellow'
rax = plt.axes([0.05, 0.7, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('40%', '60%', '80%', '100%'))
radio.on_clicked(set_hz_flag)

choose1 = plt.axes([0.3,0.05,0.1,0.075])
choose_button1 = Button(choose1,'File \n Select')
choose_button1.on_clicked(choose_transmitter)

move_to_frame = plt.axes([0.4,0.05,0.1,0.075])
move_to_frame_button = Button(move_to_frame,'mv_to_frm')
move_to_frame_button.on_clicked(move_to_frame_func)

write_csv = plt.axes([0.5,0.05,0.1,0.075])
write_csv_button = Button(write_csv,'Write \n CSV')
write_csv_button.on_clicked(write_csv_func)

runB = plt.axes([0.7,0.05,0.1,0.075])
choose_button5 = Button(runB,'RUN')
choose_button5.on_clicked(run1)


span = SpanSelector(ax3, onselect, 'horizontal',
                    rectprops=dict(alpha=0.5, facecolor='red'))


multi = MultiCursor(fig.canvas, (ax, ax2), color='r', lw=1,
                    horizOn=False, vertOn=True)


axbox_roll = plt.axes([0.07, 0.6, 0.05, 0.03])
text_box = TextBox(axbox_roll, 'Window', initial='10')
text_box.on_submit(submit_rolling)

label_box = plt.axes([0.06, 0.57, 0.15, 0.03])
label_text_box = TextBox(label_box, 'label', initial=label)
label_text_box.on_submit(change_label)

csv_box = plt.axes([0.06, 0.54, 0.15, 0.03])
csv_text_box = TextBox(csv_box, 'cvs name', initial='Output_Name')
csv_text_box.on_submit(set_csv_name)


close1 = plt.axes([0.6,0.05,0.1,0.075])
close_button = Button(close1,'close')
close_button.on_clicked(close)

plt.show()



'''
df_sub = df[['gps_lat1','gps_lat2','gps_lon1','gps_lon2','range']].dropna()

vinc = []
for bb in range(len(df_sub)):
    vinc.append(vincenty.vincenty((df_sub['gps_lat1'][bb],df_sub['gps_lon1'][bb]),(df_sub['gps_lat2'][bb],df_sub['gps_lon2'][bb])))
'''
