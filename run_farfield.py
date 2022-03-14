#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 08:39:37 2022

@author: jcerda
"""

import os
import sys
from datetime import datetime, timedelta
# from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.tamoc_plume import Plume
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_NEMO_native
from opendrift.models.ciceseoil import OpenCiceseOil
#from opendrift.models.openplume3D import OpenPlume3D

import os
import sys
import yaml
import argparse
import datetime as dt
import netCDF4 as nc
import xarray as xr
import numpy as np
from siphon.catalog import TDSCatalog
from subprocess import Popen




if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Run far field model OPENDRIFT', prog='run_farfield')
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
      endtime = (dt.datetime.today() + dt.timedelta(days=4)).strftime("%Y%m%d")

  else:
      smonth = int(PtoParam['sim']['Smonth'] + 1)
      starttime = dt.datetime(PtoParam['sim']['Syear'], PtoParam['sim']['Smonth'], PtoParam['sim']['Sday'],
                              PtoParam['sim']['Shour'], PtoParam['sim']['Smin'])

      emonth = int(PtoParam['sim']['Emonth'] + 1)
      endtime = dt.datetime(PtoParam['sim']['Eyear'], PtoParam['sim']['Emonth'], PtoParam['sim']['Eday'],
                            PtoParam['sim']['Ehour'], PtoParam['sim']['Emin'])

  if PtoParam['cicoil']['wind_factor'] != None:
      wind_factor = float(PtoParam['cicoil']['wind_factor'])
  else: wind_factor = 0.035
  
  # Select Simulation location and results output file
  output_path = output_path = os.path.join(PtoParam['outdir'],PtoParam['point']['name'])
  # TAMOC output dir
  txtTime = starttime.strftime("%Y%m%d-%H")  # str(starttime)[0:10]
  output_tamoc = os.path.join(PtoParam['outdir'],PtoParam['point']['name'],'TAMOC_output_files/')
  
  sim_name = PtoParam['sim_name']
  sim_duration = timedelta(days=int(PtoParam['cicoil']['sim_len']))
  particles_number = float(PtoParam['cicoil']['N_parti'])
  step_time = float(PtoParam['cicoil']['step_time'])    # hours
  output_step_time = float(PtoParam['cicoil']['repo_time'])    # hours
  release_points = np.int(1)

  seed_duration = endtime - starttime
  leak_days = seed_duration.days
  leak_remainder = endtime - (starttime + timedelta(days=leak_days))
  seed_days_exact = seed_duration.total_seconds() / 86400.     # seconds per day
  particles_per_day = int(np.ceil(particles_number / seed_days_exact))
  remaining_particles = int(particles_per_day * divmod(seed_days_exact, leak_days)[1])
  if remaining_particles > 0: leak_days += 1

  lon = float(PtoParam['point']['lon'])
  lat = float(PtoParam['point']['lat'])
  oil_name = PtoParam['oil']['name']
  flow_rate = float(PtoParam['spill']['flow_rate'])        # kg/s
  verticle_angle = - np.pi / 180 * int(PtoParam['spill']['vertical_angle'])   # convert degrees to radians
  horiz_angle = np.pi / 180 * int(PtoParam['spill']['horizont_angle'])      # convert degrees to radians
  GOR = float(PtoParam['spill']['GOR'])     # scf/bbl
  release_depth = float(PtoParam['spill']['depth'])       # meters
  jet_diameter = float(PtoParam['spill']['jet_diam'])         # meters
  fluid_temp = 273.15 + float(PtoParam['spill']['temp'])       # convert Celsius to Kelvin
  bins = int(PtoParam['spill']['bins'])
  
  # Create and configure OpenCiceseOil
  weathering=PtoParam['model']['weathering']
  fartype=PtoParam['model']['fartype']
  
  
  o = OpenCiceseOil(loglevel=20, weathering_model=weathering)
  if fartype == '2D Sim.':
      o.disable_vertical_motion()
  else:
      o.set_config('processes:dispersion', False)
  
  o.set_config('drift:advection_scheme',PtoParam['model']['advection_scheme'])#agregado Ang   
  o._set_config_default('drift:current_uncertainty', 0.0)
  o._set_config_default('drift:wind_uncertainty', 0.0)

  # Readers
  reader_amseas = reader_netCDF_CF_generic.Reader(filename=PtoParam['datadir'] + 'fnmoc-amseas/' + txtTime[0:8] + '/fnmoc-amseas-forecast-GoM-' + txtTime[0:8] + '-time*.nc', name='amseas_forecast')
  #print(reader_amseas)
  #reader_hycom = reader_netCDF_CF_generic.Reader(filename=PtoParam['datadir'] + 'hycom/' + 'HYCOM-forecast-GoM-' + txtTime[0:8] + '.nc', name='HYCOM_forecast') #ang
  #print(reader_hycom)
  reader_gfswinds = reader_netCDF_CF_generic.Reader(filename=PtoParam['datadir'] + 'gfs-winds/' + 'gfs-winds-forecast-GoM-' + txtTime[0:8] + '.nc', name='gfs_forecast')
  #print(reader_gfswinds)
  
  # o.add_reader([reader_basemap, reader_globcurrent, reader_oceanwind])
  o.add_reader([eval(PtoParam['model']['reader_current']),eval(PtoParam['model']['reader_winds'])]) #ang
  #o.add_reader([reader_amseas, reader_gfswinds])
  o.set_oiltype(oil_name)

  
  '''
        # Create live-oil composition
        gui.o.set_oiltype(oiltype=oil_name)
        live_comp, composition, chemdata, chemunits = gui.o.fluid_properties.get_live_composition(GOR=GOR)

        # Create vertical profiles
        for day in range(leak_days):
            tamoc_time = starttime + timedelta(days=day)
            if gui.vprofileoption.get() == 1:
                profile_name = ''.join((output_path, 'v_profile_', str(starttime)[0:10], '_day-', str(day),
                                        '_point-', str(pt), '.nc'))
                vp.create_profile(profile_name, readers[0], lon, lat, tamoc_time)

        # run near-field model
        for day in range(leak_days):
            if gui.vprofileoption.get() == 2:
                profile_name = gui.vprofilefile.var.get()
            else:
                profile_name = ''.join((output_path, 'v_profile_', str(starttime)[0:10], '_day-', str(day),
                                        '_point-', str(pt), '.nc'))
            output_name = ''.join((output_path, 'tamoc_', str(starttime)[0:10], '_day-', str(day),
                                   '_point-', str(pt)))
            gui.bpm = tamoc.run_tamoc(live_comp, composition, chemdata, chemunits, output_name, profile_name,
                                      total_flow=flow_rate, z0=release_depth, D=jet_diameter, Tj=fluid_temp,
                                      phi_0=verticle_angle, theta_0=horiz_angle, bins=bins)
  '''
  tamoc_file= ''.join((output_tamoc, 'tamoc_', txtTime, '_', PtoParam['point']['name'],'_',str(release_depth.as_integer_ratio()[0]),'m'))
  
  # Seed to Far-field and run
  spillets = particles_per_day
  for day in range(leak_days):
      tamoc_plume = Plume(plume_file=tamoc_file + '_plume.nc')
      tamoc_plume.get_particle_properties(particles_file=tamoc_file + '_particles.nc')
      
      stoptime = starttime + timedelta(days=day + 1)
      if day==leak_days - 1 and remaining_particles > 0:
          spillets = remaining_particles
          stoptime = starttime + timedelta(days=day) + leak_remainder
          
      o.seed_plume_elements(lon, lat, tamoc_plume, starttime+ timedelta(days=day), stoptime,
                               number=spillets, z_uncertainty=10, wind_drift_factor=wind_factor,
                               oiltype=oil_name)

  output_file = output_path + '/' + PtoParam['point']['name'] + '_' +sim_name + '_' + starttime.strftime("%Y%m%d-%H") + '_' + endtime.strftime("%Y%m%d-%H") + '.nc'

  o.run(duration=sim_duration,
            time_step=timedelta(hours=step_time),
            time_step_output=timedelta(hours=output_step_time),
            outfile=output_file)
  print(o)
  
  # Post processing
  postp_file=''.join((output_path,'/',PtoParam['point']['name'], '_', sim_name))
  # MAP
  o.plot(filename=postp_file +'_trayectorias.png')
  o.plot_oil_budget(filename=postp_file + '_budget.png') # agregdo por ang
  o.plot_oil_budget(show_density_viscosity=True, show_wind_and_current=True,filename=postp_file + '_budget2.png',) #ang
  # Animations
  o.animation(background=['x_sea_water_velocity', 'y_sea_water_velocity'],
              colorbar=True, fps=8, filename=postp_file + '.gif')
  o.animation(color='viscosity', fsp=8, filename=postp_file + '_viscosity.gif') #ang
  o.animation_profile(fps=8, filename=postp_file + '_profile.gif')
  
  o.animate_vertical_distribution(filename=postp_file + '_vertical.gif')
  
