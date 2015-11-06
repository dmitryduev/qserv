# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 12:13:10 2015

@author: Dr. Dmitry A. Duev
"""

import cherrypy
from cherrypy.lib import auth_digest

from jinja2 import Environment, FileSystemLoader
env = Environment(loader=FileSystemLoader('templates'))

import xml.etree.ElementTree as ET
from xml.dom import minidom
from xml.etree.ElementTree import Element

import os
import shutil

import json
from collections import OrderedDict

from dicttoxml import dicttoxml

import astropy
from astropy.io.votable import parse_single_table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def prettify(elem):
    """
        Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t")

class xmlTree(object):
    '''
        Class for handling xml files for Programs and Targets
    '''
    def __init__(self, path):
        self.path = path
    
    def getPrograms(self, programs_xml='Programs.xml'):
        self.programs_xml = programs_xml
        try:
            # parse Programs.xml:
            tree = ET.parse(os.path.join(self.path, self.programs_xml))
            self.root = tree.getroot()
        except:
            # file does not exist or empty? create 'template' then:
            with open(os.path.join(self.path, self.programs_xml), 'w') as f:
                f.write('<root>\n</root>')
            self.Programs = []
            return self.Programs
        
        # list for keeping track of programs:
        self.Programs = []

        fix_xml = False
        
        for program in self.root:
            prog = {}
            for content in program:
                prog[content.tag] = content.text
            # do not trust Programs.xml - count the number of targets!
            number_of_targets = \
                            str(self.getProgramNumberOfTargets(prog['number']))

            # build program object:
            p = Program(number=prog['number'],\
                         name=prog['name'],\
                         person_name=prog['person_name'],\
                         scientific_importance=prog['scientific_importance'],\
                         number_of_targets=number_of_targets,\
                         counter=prog['counter'],\
                         total_observation_time=prog['total_observation_time'],\
                         total_science_time=prog['total_science_time'])
            # now append:
            self.Programs.append(p)
            
            # edit Program.xml?
            if number_of_targets != prog['number_of_targets']:
                print 'Fixed wrong number of targets for program {:s}'\
                        .format(prog['name']) + ' in Programs.xml.'
                fix_xml = True

        # edit Program.xml if necessary:
        if fix_xml:
            # build new xml and write it to disk:
            tree = self.buildProgramXMLtree(self.Programs)

            with open(os.path.join(self.path, 'Programs.xml'), 'w') as f:
                f.write(minidom.parseString(ET.tostring(tree.getroot(), \
                                                'utf-8')).toprettyxml(indent="\t"))
        
        return self.Programs

    def getProgramNumberOfTargets(self, program_number):
        target_xml_path = os.path.join(self.path, \
                            'Program_{:s}'.format(program_number))
        pnot = len([f for f in os.listdir(target_xml_path) \
                        if 'Target_' in f and f[0]!='.'])
        
        return pnot
                         
    def getProgramNames(self):
        return [p.name for p in self.Programs]
        
    def getProgramNumbers(self):
        return [p.number for p in self.Programs]

    @staticmethod
    def buildProgramXMLtree(programs):
        '''
        '''
        root = ET.Element("root")
        for program in programs:
            prog = ET.SubElement(root, "Program")
            
            ET.SubElement(prog, "number").text = program.number
            ET.SubElement(prog, "name").text = program.name
            ET.SubElement(prog, "person_name").text = program.person_name
            ET.SubElement(prog, "scientific_importance").text = \
                                                program.scientific_importance
            ET.SubElement(prog, "number_of_targets").text = \
                                                program.number_of_targets
            ET.SubElement(prog, "counter").text = program.counter
            ET.SubElement(prog, "total_observation_time").text = \
                                                program.total_observation_time
            ET.SubElement(prog, "total_science_time").text = \
                                                program.total_science_time
            
        tree = ET.ElementTree(root)
        
        return tree


    def dumpProgram(self, program_number, name, number, person_name, \
                scientific_importance, number_of_targets,\
                counter, total_observation_time, total_science_time):
        '''
        '''
        try:
            self.Programs
        except:
            self.getPrograms()
        programNumbers = self.getProgramNumbers()
        
        if program_number not in programNumbers or program_number=='':
            # add new program
            self.Programs.append(Program(number, name, person_name, \
                 scientific_importance, number_of_targets, counter, \
                 total_observation_time, total_science_time))
        else:
            # change existing
            ind = [i for i,p in enumerate(self.Programs) \
                    if p.number==program_number][0]
            self.Programs[ind] = Program(number, name, person_name, \
                 scientific_importance, number_of_targets, counter, \
                 total_observation_time, total_science_time)
        
        # rename folder containing target xml files
        if program_number not in ('', number):
            os.rename(os.path.join(self.path, \
                                       'Program_{:s}'.format(program_number)),\
                      os.path.join(self.path, 'Program_{:s}'.format(number)))
        
        # build new xml and write it to disk:
        tree = self.buildProgramXMLtree(self.Programs)
        
        # make is pretty:
        with open(os.path.join(self.path, 'Programs.xml'), 'w') as f:
            f.write(minidom.parseString(ET.tostring(tree.getroot(), \
                                            'utf-8')).toprettyxml(indent="\t"))
        
        # create a directory for storing target xml-files
        target_xml_dir = os.path.join(self.path, 'Program_{:s}'.format(number))
        if not os.path.exists(target_xml_dir):
            os.makedirs(target_xml_dir)
        
        # update self.Programs:
        self.Programs = self.getPrograms()
        
        return True
        
    def removeProgram(self, name=None):
        '''
            Remove a program
        '''
        if name is None or name=='':
            return {}
        
        try:
            self.Programs
        except:
            self.getPrograms()
        
        # remove from self.Programs and Programs.xml
        try:
            program = [p for p in self.Programs if p.name==name][0]
            self.Programs.remove(program)
        except:
            return {}
        
        # build new xml and write it to disk:
        tree = self.buildProgramXMLtree(self.Programs)

        with open(os.path.join(self.path, 'Programs.xml'), 'w') as f:
            f.write(minidom.parseString(ET.tostring(tree.getroot(), \
                                            'utf-8')).toprettyxml(indent="\t"))
        
        # remove dir Program_number with target xml files
        target_xml_dir = os.path.join(self.path, \
                                        'Program_{:s}'.format(program.number))
        if os.path.exists(target_xml_dir):
            shutil.rmtree(target_xml_dir)
        
        # update self.Programs (just to make sure)
        self.Programs = self.getPrograms()
        
        return {}
        
    def removeTarget(self, program_number=None, target_number=None):
        '''
            Remove a target from a program
        '''
        if program_number is None or program_number=='' or\
            target_number is None or target_number=='':
            return {}
        
        try:
            self.Programs
        except:
            self.Programs = self.getPrograms()
        program = [p for p in self.Programs if p.number==program_number][0]
        
        target_xml = 'Target_{:s}.xml'.format(target_number)
        target_xml_path = os.path.join(self.path, \
                            'Program_{:s}'.format(program_number), target_xml)
        if os.path.exists(target_xml_path):
            os.remove(target_xml_path)
        
        # rename remaining xml files:
        for i in range(int(target_number)+1, int(program.number_of_targets)+1):
            target_xml_old = 'Target_{:d}.xml'.format(i)
            target_xml_old_path = os.path.join(self.path, \
                        'Program_{:s}'.format(program_number), target_xml_old)
            target_xml_new = 'Target_{:d}.xml'.format(i-1)
            target_xml_new_path = os.path.join(self.path, \
                        'Program_{:s}'.format(program_number), target_xml_new)
            os.rename(target_xml_old_path, target_xml_new_path)
        
        # update self.Programs
        self.Programs = self.getPrograms()
        
        return {}
        
        
    def getTargets(self, program, target_list_xml=None):
        '''
            Get a list of targets (each being a dict) for a given program
        '''
        if target_list_xml is None:
            target_list_xml = ['Target_{:d}.xml'.format(i+1) \
                                for i in range(int(program.number_of_targets))]
        
        # list for keeping track of targets for each program:
        try:
            self.Targets
        except:
            self.Targets = {}
        
        self.Targets[program.name] = []
        
        for target_xml in target_list_xml:
            tree = ET.parse(os.path.join(self.path, \
                            'Program_{:s}'.format(program.number), target_xml))
            root = tree.getroot()
        
            targ = {}
            targ['Object'] = []
            for content in root:
                if content.tag != 'Object':
                    targ[content.tag] = content.text
                else:
                    obj = {}
                    obj['Observation'] = []
                    for data_obj in content:
                        if data_obj.tag != 'Observation':
                            obj[data_obj.tag] = data_obj.text
                        else:
                            obs = {}
                            for data_obs in data_obj:
                                obs[data_obs.tag] = data_obs.text
                            obj['Observation'].append(obs)
                    targ['Object'].append(obj)

            self.Targets[program.name].append(targ)
#        print 'N_targ = {:d}'.format(len(self.Targets[program.name]))
        return self.Targets[program.name]


    def getTargetNames(self, program):
        '''
            Get target names for program
            
            Each target must have a unique name!
        '''
        try:
            self.Targets[program.name]
        except:
            self.getTargets(program)
        
        targetNames = [t['name'] for t in self.Targets[program.name]]
        
        return targetNames
        
        
    def batchEditTargets(self, program, time_critical_flag="", 
                          visited_times_for_completion="", seeing_limit="",
                          cadence="", obj_epoch="", 
                          obj_sun_altitude_limit="",
                          obj_moon_phase_window="",
                          obj_airmass_limit="",
                          obj_sun_distance_limit="",
                          obj_moon_distance_limit="",
                          obj_sky_brightness_limit="",
                          obj_hour_angle_limit="",
                          obs_exposure_time="", obs_ao_flag="", 
                          obs_filter_code="", obs_repeat_times=""):
        target_list_xml = ['Target_{:d}.xml'.format(i+1) \
                            for i in range(int(program.number_of_targets))]
        
        for target_xml in target_list_xml:
            target_xml_path = os.path.join(self.path, \
                            'Program_{:s}'.format(program.number), target_xml)
            tree = ET.parse(target_xml_path)
            root = tree.getroot()

            if time_critical_flag!="":
                # does the tag exist? if not, create
                if root.find('time_critical_flag') is None:
                    root.append(Element('time_critical_flag'))
                tag = root.find('time_critical_flag')
                tag.text = time_critical_flag
                
            if visited_times_for_completion!="":
                # does the tag exist? if not, create
                if root.find('visited_times_for_completion') is None:
                    root.append(Element('visited_times_for_completion'))
                tag = root.find('visited_times_for_completion')
                tag.text = visited_times_for_completion
                
            if seeing_limit!="":
                # does the tag exist? if not, create
                if root.find('seeing_limit') is None:
                    root.append(Element('seeing_limit'))
                tag = root.find('seeing_limit')
                tag.text = seeing_limit
                
            if cadence!="":
                # does the tag exist? if not, create
                if root.find('cadence') is None:
                    root.append(Element('cadence'))
                tag = root.find('cadence')
                tag.text = cadence
            
            # iterate over Objects:
            objs = root.findall('Object')
            for obj in objs:
                if obj_epoch!="":
                    # does the tag exist? if not, create
                    if obj.find('epoch') is None:
                        obj.append(Element('epoch'))
                    tag = obj.find('epoch')
                    tag.text = obj_epoch

                if obj_sun_altitude_limit!="":
                    # does the tag exist? if not, create
                    if obj.find('sun_altitude_limit') is None:
                        obj.append(Element('sun_altitude_limit'))
                    tag = obj.find('sun_altitude_limit')
                    tag.text = obj_sun_altitude_limit
                    
                if obj_moon_phase_window!="":
                    # does the tag exist? if not, create
                    if obj.find('moon_phase_window') is None:
                        obj.append(Element('moon_phase_window'))
                    tag = obj.find('moon_phase_window')
                    tag.text = obj_moon_phase_window
                    
                if obj_airmass_limit!="":
                    # does the tag exist? if not, create
                    if obj.find('airmass_limit') is None:
                        obj.append(Element('airmass_limit'))
                    tag = obj.find('airmass_limit')
                    tag.text = obj_airmass_limit
                    
                if obj_sun_distance_limit!="":
                    # does the tag exist? if not, create
                    if obj.find('sun_distance_limit') is None:
                        obj.append(Element('sun_distance_limit'))
                    tag = obj.find('sun_distance_limit')
                    tag.text = obj_sun_distance_limit
                    
                if obj_moon_distance_limit!="":
                    # does the tag exist? if not, create
                    if obj.find('moon_distance_limit') is None:
                        obj.append(Element('moon_distance_limit'))
                    tag = obj.find('moon_distance_limit')
                    tag.text = obj_moon_distance_limit
                    
                if obj_sky_brightness_limit!="":
                    # does the tag exist? if not, create
                    if obj.find('sky_brightness_limit') is None:
                        obj.append(Element('sky_brightness_limit'))
                    tag = obj.find('sky_brightness_limit')
                    tag.text = obj_sky_brightness_limit
                    
                if obj_hour_angle_limit!="":
                    # does the tag exist? if not, create
                    if obj.find('hour_angle_limit') is None:
                        obj.append(Element('hour_angle_limit'))
                    tag = obj.find('hour_angle_limit')
                    tag.text = obj_hour_angle_limit
                    
                # iterate over Observations:
                obss = obj.findall('Observation')
                for obs in obss:
                    if obs_exposure_time!="":
                        # does the tag exist? if not, create
                        if obs.find('exposure_time') is None:
                            obs.append(Element('exposure_time'))
                        tag = obs.find('exposure_time')
                        tag.text = obs_exposure_time

                    if obs_ao_flag!="":
                        # does the tag exist? if not, create
                        if obs.find('ao_flag') is None:
                            obs.append(Element('ao_flag'))
                        tag = obs.find('ao_flag')
                        tag.text = obs_ao_flag
                
                    if obs_filter_code!="":
                        # does the tag exist? if not, create
                        if obs.find('filter_code') is None:
                            obs.append(Element('filter_code'))
                        tag = obs.find('filter_code')
                        tag.text = obs_filter_code
                
                    if obs_repeat_times!="":
                        # does the tag exist? if not, create
                        if obs.find('repeat_times') is None:
                            obs.append(Element('repeat_times'))
                        tag = obs.find('repeat_times')
                        tag.text = obs_repeat_times
                

            # save updated file:
#            with open(target_xml_path, 'w') as f:
#                out = minidom.parseString(ET.tostring(tree.getroot(), \
#                    'utf-8').replace(' ', ' ')).toprettyxml(indent='', newl='')
#                f.write(out.replace('<?xml version="1.0" ?>', ''))

            # build an xml-file:
            # this is good enough, but adds unnecessary <item> tags. remove em:
            target_xml = minidom.parseString(ET.tostring(tree.getroot(), \
                                             'utf-8')).toprettyxml()
            # <item>'s left extra \t's after them - remove them:
#            target_xml = target_xml.replace('\t\t\t','\t\t')
#            target_xml = target_xml.replace('\t\t\t\t','\t\t\t')
            target_xml = target_xml.replace('<?xml version="1.0" ?>', '')
            target_xml = target_xml.split('\n')
            target_xml = [t for t in target_xml if 'item>' not in t and
                            not t.isspace()]


#            print target_xml_path
            with open(target_xml_path, 'w') as f:
                for line in target_xml[1:-1]:
                    f.write('{:s}\n'.format(line))
                f.write('{:s}'.format(target_xml[-1]))
    
    def dumpTarget(self, program, target_number, target):
        '''
            Edit or create target xml file
        '''
        if target_number=='':
            targets_added = 1
            # program not empty?
            if int(program.number_of_targets)!=0:
                # get targets
                targets = self.getTargets(program)
                
                # get max target number. when adding new tragets, start from max_number+1
                max_number = max([int(t['number']) for t in targets])
            else:
                max_number = 0
            target_number = max_number+1
            target['number'] = str(target_number)
        else:
            targets_added = 0
            target_number = int(target_number)
        # file Target_*.xml number need not be = target_number from the xml file
        xml_file_number = target_number
            
        # build an xml-file:

        target_xml = dicttoxml(target, custom_root='Target', attr_type=False)
        # this is good enough, but adds unnecessary <item> tags. remove em:
        dom = minidom.parseString(target_xml)
        target_xml = dom.toprettyxml()
        # <item>'s left extra \t's after them - remove them:
        target_xml = target_xml.replace('\t\t\t','\t\t')
        target_xml = target_xml.replace('\t\t\t\t','\t\t\t')
        target_xml = target_xml.replace('<?xml version="1.0" ?>', '')
        target_xml = target_xml.split('\n')
        target_xml = [t for t in target_xml if 'item>' not in t]
        # deal with missing <Object>s and <Observation>s:
#        xml_out = []
#        for line in target_xml[1:-1]:
#            xml_out.append('{:s}\n'.format(line))
#        xml_out.append('{:s}'.format(target_xml[-1]))
        ind_obs_start = [i for i,v in enumerate(target_xml) if '<Observation>' in v]
        ind_obs_stop = [i for i,v in enumerate(target_xml) if '</Observation>' in v]
        for (start, stop) in zip(ind_obs_start, ind_obs_stop):
            ind_num_obs = [i+start for i,v in enumerate(target_xml[start:stop]) \
                                if '<number>' in v]
            if len(ind_num_obs)>1:
                for ind in ind_num_obs[:0:-1]:
                    target_xml.insert(ind, '\t\t</Observation>\n\t\t<Observation>')

        ind_obj= [i for i,v in enumerate(target_xml) if v[:10]=='\t\t<number>']
        for ind in ind_obj[:0:-1]:
            target_xml.insert(ind, '\t</Object>\n\t<Object>')

#        print target_xml
        
        target_xml_path = os.path.join(self.path, \
                        'Program_{:s}'.format(program.number), 
                        'Target_{:d}.xml'.format(xml_file_number))

        with open(target_xml_path, 'w') as f:
            for line in target_xml[1:-1]:
                f.write('{:s}\n'.format(line))
            f.write('{:s}'.format(target_xml[-1]))

        # update program number of targets!!
        self.dumpProgram(program.number, \
                program.name, program.number, program.person_name, \
                program.scientific_importance, \
                str(int(program.number_of_targets)+targets_added),\
                program.counter, program.total_observation_time, \
                program.total_science_time)
                
        cherrypy.log('Target_{:d}.xml edited/created in Program_{:s}'.\
                format(xml_file_number, program.number))
        
    
    def dumpTargetList(self, program, data):
        '''
            List from an external file
        '''
        # program not empty?
        if int(program.number_of_targets)!=0:
            # get targets
            targets = self.getTargets(program)
            
            # get max target number. when adding new tragets, start from max_number+1
            max_number = max([int(t['number']) for t in targets])
        else:
            max_number = 0
        
        # get existing target names:
        targetNames = self.getTargetNames(program)
        
        targets_added = 0
        
        # do some guessing about the table
        if type(data) == astropy.io.votable.tree.Table:
            # VOtable? convert to normal table
            table = data.to_table()
        else:
            table = data
        # target names:
        target_name_field = [table[f].name for f in table.colnames \
                                if 'name' in table[f].name or \
                                (table[f].description is not None and\
                                'name' in table[f].description)][0]
        if type(table[target_name_field].data) is not np.ndarray:
            target_names_list = table[target_name_field].data.data
        else:
            target_names_list = table[target_name_field].data
        # RA/Dec:
        if '_RAJ2000' in table.colnames and '_DEJ2000' in table.colnames:
            # no columns? add them:
            ra = [ra.replace(' ', ':') for ra in table['_RAJ2000']]
            dec = [dec.replace(' ', ':') for dec in table['_DEJ2000']]
        else:
            # ra/dec in degrees?
            if not 'RAJ2000' in table.colnames or not 'DEJ2000' in table.colnames:
                raise Exception('Could not find coordinates in the imported file')
            else:
                # get units:
                if table['RAJ2000'].description is not None and \
                        not 'degree' in table['RAJ2000'].description:
                    raise Exception('RAJ2000 must be in degrees!')
                else:
                    crd = SkyCoord(ra=ra, dec=dec, \
                        unit=(u.deg, u.deg), frame='icrs').to_string('hmsdms')
                    crd_str = np.array([c.split() for c in crd])
                    ra = [ra.replace('h',':').replace('m',':').replace('s','')\
                            for ra in crd_str[:,0]]
                    dec = [dec.replace('d',':').replace('m',':').replace('s','')\
                            for dec in crd_str[:,1]]
        # no comments?
        if 'comment' in table.colnames:
            if type(table['comment'].data) is not np.ndarray:
                comment = table['comment'].data.data
            else:
                comment = table['comment'].data
        else:
            comment = ['']*len(target_names_list)
        # Vmag/Vemag/mag
        if 'Vmag' in table.colnames:
            mag = table['Vmag']
        elif 'Vemag' in table.colnames:
            mag = table['Vemag']
        elif 'mag' in table.colnames:
            mag = table['mag']
            
        # epoch should be J2000?
        if 'epoch' in table.colnames:
            epoch = table['epoch']
        else:
            epoch = ['2000.0']*len(target_names_list)
        
        # iterate over entries in the parsed table
        for ei, _ in enumerate(table):            
            # target name must be unique! check it and skip entry if necessary:
            if target_names_list[ei] in targetNames:
                continue
            # Ordnung muss sein!
            target = OrderedDict([('program_number',program.number),
                      ('number',str(max_number+ei+1)),
                      ('name',str(target_names_list[ei])), \
                      ('time_critical_flag','0'),
                      ('visited_times_for_completion','1'),
                      ('seeing_limit',''), ('visited_times','0'), ('done','0'),
                      ('cadence','0'), ('comment',str(comment[ei])),
                      ('Object',[])])
            target['Object'].append(OrderedDict([('number','1'), 
                                 ('RA',ra[ei]), ('dec',dec[ei]),
                                 ('epoch',epoch[ei]), ('magnitude',mag[ei]),
                                 ('sun_altitude_limit',''), ('moon_phase_window',''),
                                 ('airmass_limit',''), ('sun_distance_limit',''),
                                 ('moon_distance_limit',''), ('sky_brightness_limit',''),
                                 ('hour_angle_limit',''), ('done','0'),
                                 ('Observation',[])]))
            for oi, obj in enumerate(target['Object']):
                target['Object'][oi]['Observation'].append(\
                        OrderedDict([('number','1'),
                        ('exposure_time','0'),('ao_flag','1'), ('filter_code',''),
                        ('camera_mode',''),('repeat_times','0'),
                        ('repeated','0'),('done','0')]))

            # build an xml-file:
            target_xml = dicttoxml(target, custom_root='Target', attr_type=False)
            # this is good enough, but adds unnecessary <item> tags. remove em:
            dom = minidom.parseString(target_xml)
            target_xml = dom.toprettyxml()
            # <item>'s left extra \t's after them - remove them:
            target_xml = target_xml.replace('\t\t\t','\t\t')
            target_xml = target_xml.replace('\t\t\t\t','\t\t\t')
            target_xml = target_xml.replace('<?xml version="1.0" ?>', '')
            target_xml = target_xml.split('\n')
            target_xml = [t for t in target_xml if 'item>' not in t]
            
#            print target_xml
            target_xml_path = os.path.join(self.path, \
                            'Program_{:s}'.format(program.number), 
                            'Target_{:d}.xml'.format(\
                                    int(program.number_of_targets)+ei+1))
#            print target_xml_path
            with open(target_xml_path, 'w') as f:
                for line in target_xml[1:-1]:
                    f.write('{:s}\n'.format(line))
                f.write('{:s}'.format(target_xml[-1]))
                
            # count added targets:
            targets_added += 1
                
        # update program number of targets!!
        self.dumpProgram(program.number, \
                program.name, program.number, program.person_name, \
                program.scientific_importance, \
                str(int(program.number_of_targets)+targets_added),\
                program.counter, program.total_observation_time, \
                program.total_science_time)

#@cherrypy.popargs('name')
class Root(object):
    
    def __init__(self, path_to_queue):
        self.path_to_queue = path_to_queue
    
    # URL dispatcher
    def _cp_dispatch(self, vpath):
        if len(vpath) == 1:
            cherrypy.request.params['prog_number'] = vpath.pop()
            return self

#        if len(vpath) == 2:
#            cherrypy.request.params['targets'] = vpath.pop(0)  # /band name/
#            return self

        return vpath

    
    @cherrypy.expose
    def index(self, prog_number=None):
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # get entries:
        programs = xmlT.getPrograms(programs_xml='Programs.xml')
        
        if prog_number is None:
            # render programs:
            tmpl = env.get_template('index.html')
            return tmpl.render(programs=programs)
        # display Programs.xml:
        elif prog_number == 'Programs.xml':
            # render Programs.xml:
            with open(os.path.join(xmlT.path, prog_number), 'r') as f:
                page = ''.join('{:s}'.format(line) for line in f)
                page = '<pre>'+page+'</pre>'
            return page
#        # display Target_n.xml:
#        elif 'Program' in prog_number and 'Target' in prog_number \
#                and 'xml' in prog_number:
#            pass
#            page = ''
#            return page
        else:
            #render targets
            # get program:
            try:
                number = prog_number.split('_')[1]
                program = [p for p in programs if p.number==number][0]
            except:
                raise Exception('Program_{:s} not found.'.format(number))
#            from time import time as _time
#            tic = _time()
            targets = xmlT.getTargets(program)
#            print 'getting targets took {:f} seconds'.format(_time()-tic)
            # populate template
#            tic = _time()
            tmpl = env.get_template('targets.html')
#            print 'loading template took {:f} seconds'.format(_time()-tic)
#            tic = _time()
#            tmpl.render(targets=targets, programName=program.name,
#                               programNumber=program.number)
#            print 'rendering template took {:f} seconds'.format(_time()-tic)
            return tmpl.render(targets=targets, programName=program.name,
                               programNumber=program.number)
    
    # request and receive Program params in json format
    @cherrypy.expose
    def prog_param(self, program_name=None):
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # get entries:
        programs = xmlT.getPrograms(programs_xml='Programs.xml')
        
        if program_name is not None:
            # construct json object:
            prog = [p for p in programs if p.name==program_name]
            # found?
            if len(prog)>0:
                json_obj = prog[0].makeJSON()
                return json_obj
            else:
                return {}
        else:
            return {}
            
    # request and receive Target params in json format
    @cherrypy.expose
    def targ_param(self, program_number=None, target_number=None):
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # get entries:
        programs = xmlT.getPrograms(programs_xml='Programs.xml')
        
        if program_number is not None and \
                (target_number is not None and target_number!=""):
            # construct json object:
            prog = [p for p in programs if p.number==program_number]
            # found?
            if len(prog)>0:
                target_xml = ['Target_{:d}.xml'.format(int(target_number))]
                target_dict = xmlT.getTargets(program=prog[0], \
                                                target_list_xml=target_xml)[0]
                json_obj = json.dumps(target_dict)
#                print json_obj
                return json_obj
            else:
                return {}
        else:
            return {}
                
    # save new/edited Program
    @cherrypy.expose
    def save(self, program_number=None, \
                name=None, number=None, person_name=None, \
                scientific_importance=None, number_of_targets=None,\
                counter=None, total_observation_time=None, \
                total_science_time=None):
        # bad input:
        if None in (program_number, name, number, person_name, \
                scientific_importance, number_of_targets,\
                counter, total_observation_time, total_science_time) or\
            (name=='' or number=='' or person_name=='' or\
             scientific_importance==''):
            return {}
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # read in Programs:
        xmlT.getPrograms(programs_xml='Programs.xml')
        # save program:
        xmlT.dumpProgram(program_number, name, number, person_name, \
                scientific_importance, number_of_targets,\
                counter, total_observation_time, total_science_time)

        return {}
        
    # remove new/edited Program
    @cherrypy.expose
    def remove(self, program_name=None):
        # bad input:
        if program_name is None:
            return {}
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # read in Programs:
        xmlT.getPrograms(programs_xml='Programs.xml')
        # get program names:
        xmlT.removeProgram(program_name)

        return {}
    
    @cherrypy.expose
    def targetBatchUpdate(self, program_number, time_critical_flag="", 
                          visited_times_for_completion="", seeing_limit="",
                          cadence="", obj_epoch="", 
                          obj_sun_altitude_limit="",
                          obj_moon_phase_window="",
                          obj_airmass_limit="",
                          obj_sun_distance_limit="",
                          obj_moon_distance_limit="",
                          obj_sky_brightness_limit="",
                          obj_hour_angle_limit="",
                          obs_exposure_time="", obs_ao_flag="", 
                          obs_filter_code="", obs_repeat_times=""):
        
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # get entries:
        programs = xmlT.getPrograms(programs_xml='Programs.xml')
        
        program = [p for p in programs if p.number==program_number][0]
        
        xmlT.batchEditTargets(program, time_critical_flag, 
                          visited_times_for_completion, seeing_limit,
                          cadence, obj_epoch, obj_sun_altitude_limit,
                          obj_moon_phase_window, obj_airmass_limit,
                          obj_sun_distance_limit, obj_moon_distance_limit,
                          obj_sky_brightness_limit, obj_hour_angle_limit,
                          obs_exposure_time, obs_ao_flag, 
                          obs_filter_code, obs_repeat_times)
        
        return {}
        
    @cherrypy.expose
    def importTargetList(self, targetList, program_number):
        # read in in chunks:
#        size = 0
#        while True:
#            data = targetList.file.read(8192)
#            if not data:
#                break
#            size += len(data)

#        # read everything in one go:
#        data = targetList.file.readlines()
#        # skip comments:
#        data = [d for d in data if d[0]!='#']

        # let's build an astropy table following Vizier's
        # column naming convention
        
        # is it a VOtable?
        if targetList.filename[-3:]=='vot':
            data = parse_single_table(targetList.file)
        # if not - it must be readable by astropy.io.ascii
        else:
            data = ascii.read(targetList.file)

        # load Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # get entries:
        programs = xmlT.getPrograms(programs_xml='Programs.xml')
        
        program = [p for p in programs if p.number==program_number][0]
        
        xmlT.dumpTargetList(program, data)
        
        raise cherrypy.HTTPRedirect('/Program_{:s}'.format(program_number)) 
        
    @cherrypy.expose
    def targetUpdate(self, **kargs):
        '''
            Update Target xml file
        '''
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # get entries:
        programs = xmlT.getPrograms(programs_xml='Programs.xml')
        
        program = [p for p in programs if p.number==kargs['program_number']][0]
        
        
        nObj = len([k for k,v in kargs.iteritems() if 'obj_number' in k])
        nObs = [len([k for k,v in kargs.iteritems() \
                    if 'obs_number_{:d}'.format(i+1) in k]) for i in range(nObj)]
#        print nObj,nObs
        # make target dict:
        target = OrderedDict((('program_number',kargs['program_number']),\
                  ('number',kargs['number']),\
                  ('name',kargs['name']),\
                  ('visited_times_for_completion',kargs['visited_times_for_completion']),\
                  ('seeing_limit',kargs['seeing_limit']),\
                  ('visited_times',kargs['visited_times']),\
                  ('done',kargs['done']),\
                  ('cadence',kargs['cadence']),\
                  ('comment',kargs['comment']),\
                  ('time_critical_flag',kargs['time_critical_flag']),\
                  ('Object',[])))
        obj_numbers = sorted([s[-s[::-1].index('_'):] for s in kargs if 'obj_number' in s])
        obs_numbers = [ sorted([s[-s[::-1].index('_'):] for s in kargs \
                        if 'obs_number_{:d}'.format(ii+1) in s]) for ii in range(nObj)]

        for nOj in range(nObj):
            # fix number if necessary
            target['Object'].append(OrderedDict(( \
                ('number',nOj+1),\
                ('RA',kargs['obj_RA_{:s}'.format(obj_numbers[nOj])]),\
                ('dec',kargs['obj_dec_{:s}'.format(obj_numbers[nOj])]),\
                ('ra_rate',kargs['obj_ra_rate_{:s}'.format(obj_numbers[nOj])]),\
                ('dec_rate',kargs['obj_dec_rate_{:s}'.format(obj_numbers[nOj])]),\
                ('epoch',kargs['obj_epoch_{:s}'.format(obj_numbers[nOj])]),\
                ('magnitude',kargs['obj_magnitude_{:s}'.format(obj_numbers[nOj])]),\
                ('sun_altitude_limit',\
                    kargs['obj_sun_altitude_limit_{:s}'.format(obj_numbers[nOj])]),\
                ('moon_phase_window',\
                    kargs['obj_moon_phase_window_{:s}'.format(obj_numbers[nOj])]),\
                ('airmass_limit',\
                    kargs['obj_airmass_limit_{:s}'.format(obj_numbers[nOj])]),\
                ('sun_distance_limit',\
                    kargs['obj_sun_distance_limit_{:s}'.format(obj_numbers[nOj])]),\
                ('moon_distance_limit',\
                    kargs['obj_moon_distance_limit_{:s}'.format(obj_numbers[nOj])]),\
                ('sky_brightness_limit',\
                    kargs['obj_sky_brightness_limit_{:s}'.format(obj_numbers[nOj])]),\
                ('hour_angle_limit',\
                    kargs['obj_hour_angle_limit_{:s}'.format(obj_numbers[nOj])]),\
                ('done',kargs['obj_done_{:s}'.format(obj_numbers[nOj])]),\
                ('Observation',[])
            )))
            for nOs in range(nObs[nOj]):
                # fix number if necessary
                target['Object'][nOj]['Observation'].append(OrderedDict(( \
                    ('number',nOs+1),\
                    ('exposure_time',kargs['obs_exposure_time_{:s}_{:s}'.\
                            format(obj_numbers[nOj], obs_numbers[nOj][nOs])]),\
                    ('ao_flag',kargs['obs_ao_flag_{:s}_{:s}'.\
                            format(obj_numbers[nOj], obs_numbers[nOj][nOs])]),\
                    ('filter_code',kargs['obs_filter_code_{:s}_{:s}'.\
                            format(obj_numbers[nOj], obs_numbers[nOj][nOs])]),\
                    ('camera_mode',kargs['obs_camera_mode_{:s}_{:s}'.\
                            format(obj_numbers[nOj], obs_numbers[nOj][nOs])]),\
                    ('repeat_times',kargs['obs_repeat_times_{:s}_{:s}'.\
                            format(obj_numbers[nOj], obs_numbers[nOj][nOs])]),\
                    ('repeated',kargs['obs_repeated_{:s}_{:s}'.\
                            format(obj_numbers[nOj], obs_numbers[nOj][nOs])]),\
                    ('done',kargs['obs_done_{:s}_{:s}'.\
                            format(obj_numbers[nOj], obs_numbers[nOj][nOs])])
                )))

#        print target
        
        xmlT.dumpTarget(program, kargs['target_number'], target)
        
        # save new/edited Program
    @cherrypy.expose
    def removeTarget(self, program_number=None, target_number=None):
        # bad input:
        if program_number is None or program_number=="" \
            or target_number is None or target_number=="":
            return {}
        # read in Programs.xml:
        path = self.path_to_queue
        xmlT = xmlTree(path)
        # read in Programs:
        xmlT.getPrograms(programs_xml='Programs.xml')
        # get program names:
        xmlT.removeTarget(program_number, target_number)
        
        cherrypy.log('removed Target_{:s}.xml from Program_{:s}'.\
                        format(target_number, program_number))
        cherrypy.log('Note that remaining target xml files were ranaimed if '+\
                'target_number<number_of_targets to keep file numbering order')
        
        return {}
        

class Program(object):
    def __init__(self, number, name, person_name, scientific_importance, \
                 number_of_targets, counter, total_observation_time, \
                 total_science_time):
        self.number = number
        self.name = name
        self.person_name = person_name
        self.scientific_importance = scientific_importance
        self.number_of_targets = number_of_targets
        self.counter = counter
        self.total_observation_time = total_observation_time
        self.total_science_time = total_science_time

    def makeJSON(self):
        dic = OrderedDict([['number', self.number],\
                           ['name', self.name],\
                           ['person_name', self.person_name],\
                           ['scientific_importance', self.scientific_importance],\
                           ['number_of_targets', self.number_of_targets],\
                           ['counter', self.counter],\
                           ['total_observation_time', self.total_observation_time],\
                           ['total_science_time', self.total_science_time]])
        return json.dumps(dic)



if __name__ == '__main__':
#    cherrypy.quickstart(Root())
    
    USERS = {'admin': 'robo@0'}
    
    cherrypy.config.update({'server.socket_host': '0.0.0.0',
                             'server.socket_port': 8080,
                             'server.thread_pool' : 8,
                             'log.access_file': 'server_access.log',
                             'log.error_file': 'server_actions.log'
                            })

    conf = {
         '/': {
             'tools.sessions.on': True,
             'tools.staticdir.root': os.path.abspath(os.getcwd()),
             'tools.auth_digest.on': True,
             'tools.auth_digest.realm': 'hola!',
             'tools.auth_digest.get_ha1': auth_digest.get_ha1_dict_plain(USERS),
             'tools.auth_digest.key': 'd8765asdf6c787ag333'
         },
         '/static': {
             'tools.staticdir.on': True,
#             'tools.staticdir.dir': os.path.join(os.path.abspath(os.getcwd()), 'public')
             'tools.staticdir.dir': './public',
             'tools.auth_digest.on': True,
             'tools.auth_digest.realm': 'hola!',
             'tools.auth_digest.get_ha1': auth_digest.get_ha1_dict_plain(USERS),
             'tools.auth_digest.key': 'd8765asdf6c787ag333'
         }
    }
#    path_to_queue = './'
    path_to_queue = '/Users/dmitryduev/_caltech/roboao/Queue/'
#    path_to_queue = '/Users/dmitryduev/web/qserv/test/'
    cherrypy.quickstart(Root(path_to_queue), '/', conf)