#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 09:59:53 2022

@author: jcerda
Tomado: Koztaz GUI
"""
import os
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.ciceseoil import OpenCiceseOil
import opendrift.models.vertical_profile as vp
import opendrift.models.run_tamoc as tamoc
import numpy as np
import yaml
import argparse
import datetime as dt


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Run near field model TAMOC', prog='run_nearfield')
  parser.add_argument('--subset', '-s', action='store', dest='PtoParam', help='yaml file with the information about point.')
  parser.add_argument('--commands', '-c', action='store_true', dest='show_commands', help='Just show the commands to run')

  # args = parser.parse_args('--config-file', "../Pto_Config.yaml")
  args = parser.parse_args()
  # TODO Add verbose mode
  if args.PtoParam:
      with open(args.PtoParam, 'r') as stream:
          try: 
              PtoParam = yaml.safe_load(stream)
          except:
              print ('Something went wrong reading ' + args.subsetconfig)

  # IF date is NOT given, run the present day
  if PtoParam['sim']['StarTime'] == None:
      starttime = (dt.datetime.today() - dt.timedelta(days=1)).strftime("%Y%m%d")
  else:
      starttime = dt.datetime(PtoParam['sim']['Syear'], PtoParam['sim']['Smonth'], PtoParam['sim']['Sday'],
                           PtoParam['sim']['Shour'], PtoParam['sim']['Smin'])

  # output paths
  txtTime = starttime.strftime("%Y%m%d-%H")  # str(starttime)[0:10]
  output_path = os.path.join(PtoParam['outdir'],PtoParam['point']['name'],'TAMOC_output_files/')
  try:
      os.mkdir(output_path) 
  except OSError as error:
      print(error)
  
  # load Point info
  lon = float(PtoParam['point']['lon'])
  lat = float(PtoParam['point']['lat'])
  
  # load oil properties
  oil_origin = PtoParam['oil']['origin'] # Ex: "GULF OF MEXICO, USA"
  oil_name = PtoParam['oil']['name']     # Ex: "GENERIC BUNKER C"
  flow_rate = float(PtoParam['spill']['flow_rate'])  # Ex: 40.0 kg/s
  verticle_angle = - np.pi / 180 * float(PtoParam['spill']['vertical_angle'])  # Convert degrees to radians
  horiz_angle = np.pi / 180 * float(PtoParam['spill']['horizont_angle'])       # Convert degrees to radians
  GOR = float(PtoParam['spill']['GOR'])  # Ex: 1500.0 scf/bbl
  release_depth = float(PtoParam['spill']['depth'])  # Ex: 400.0 m
  jet_diameter = float(PtoParam['spill']['jet_diam'])  # Ex: 0.2 m
  fluid_temp = 273.15 + float(PtoParam['spill']['temp'])  # Ex: 65.0 convert Celsius to Kelvin
  bins = int(PtoParam['spill']['bins'])  # Ex: 12
  
   # Create live-oil composition
  weathering=PtoParam['model']['weathering']
  llevel=PtoParam['model']['log_level']
  o = OpenCiceseOil(loglevel=llevel, weathering_model=weathering)
  o.set_oiltype(oiltype=oil_name)
  live_comp, composition, chemdata, chemunits = o.fluid_properties.get_live_composition(GOR=GOR)

  # Create vertical profiles
  tamoc_time = starttime
  ncDataIn= ''.join([PtoParam['datadir'],'fnmoc-amseas/',txtTime[0:8],'/fnmoc-amseas-forecast-GoM-', txtTime[0:8], '-time',txtTime[9:],'.nc'])
  profile_name = ''.join((output_path, 'v_profile_', txtTime, '_', PtoParam['point']['name'], '.nc'))

  readers=reader_netCDF_CF_generic.Reader(ncDataIn)
  vp.create_profile(profile_name, readers, lon, lat, tamoc_time)

  # run near-field model
  output_name = ''.join((output_path, 'tamoc_', txtTime, '_', PtoParam['point']['name'],'_',str(release_depth.as_integer_ratio()[0]),'m'))
  tamoc.run_tamoc(live_comp, composition, chemdata, chemunits, output_name, profile_name,
                                  total_flow=flow_rate, z0=release_depth, D=jet_diameter, Tj=fluid_temp,
                                  phi_0=verticle_angle, theta_0=horiz_angle, bins=bins)
