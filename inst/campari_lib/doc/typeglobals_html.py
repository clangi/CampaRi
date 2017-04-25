#--------------------------------------------------------------------------#
# LICENSE INFO:                                                            #
#--------------------------------------------------------------------------#
#    This file is part of CAMPARI.                                         #
#                                                                          #
#    Version 2.0                                                           #
#                                                                          #
#    Copyright (C) 2014, The CAMPARI development team (current and former  #
#                        contributors)                                     #
#                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang #
#                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     #
#                        Nicholas Lyle, Nicolas Bloechliger                #
#                                                                          #
#    Website: http://sourceforge.net/projects/campari/                     #
#                                                                          #
#    CAMPARI is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by  #
#    the Free Software Foundation, either version 3 of the License, or     #
#    (at your option) any later version.                                   #
#                                                                          #
#    CAMPARI is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#    GNU General Public License for more details.                          #
#                                                                          #
#    You should have received a copy of the GNU General Public License     #
#    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      #
#--------------------------------------------------------------------------#
# AUTHORSHIP INFO:                                                         #
#--------------------------------------------------------------------------#
#                                                                          #
# MAIN AUTHOR:   Adam Steffen                                              #
#                                                                          #
#--------------------------------------------------------------------------#

import string
import re
import os

#################################################################################
#
# Builds a nested set of hash tables for all mod_* files in the current directory
# The organization of the tables is:
#
#      1.  File Name (mod_*)
#         A. Variable Type (integer,logical,...)
#            1) a variable of that time
#            2) another variable of that type
#           ...
#            n)
#         B.
#        ...
#      2.
#     ...
#################################################################################

# build array of all globals and store routines that have access, access and modify the global.
globals = {}

# Files (mod_*)
modules = {}
descriptions = {}
tescriptions = {}
customtypes = {}
listofcustomtypes = {}

# loop over every mod_* file in current directory to extract custom types
for file in os.listdir("."):
  if (file.startswith('mod_') and not file.endswith('_interfaces.f90')):
    # open each mod_* file
    fileIN = open(file)
    line = fileIN.readline().lstrip(' ').rstrip()
    # Types (integer,etc..)
    module = {}
    # list of types temporary
    listoftypes = {}
    flag1 = False # add only variables inside module block
    # loop over all lines in module block
    while line:
      # build hash table of global descriptions
      if re.match("!.*?\:",line):
        type, descpt = line.split(':',1)
        type = type + " "
        x, name, y = re.split('\s+',type,2)
        descriptions[name] = descpt
      # reach end of module block
      if line.startswith("end module"):
        break
      # start adding variables
      elif line.startswith("module "):
        flag1 = True
      # add all variables from a type block
      # temporary new Type is created and removed at the end of the block
      elif line.startswith("type ") and not line.startswith("type (") and not line.startswith("type  ("):
        x, metatype = re.split('\s+',line,1)
        listofcustomtypes[metatype] = metatype
        descriptor = ""
        line = fileIN.readline()
        if line.count(" !") > 0:
          line, descriptor = line.lstrip().rstrip().split(" !",1)
          descriptor = descriptor.lstrip().rstrip()
        else:
          line = line.lstrip().rstrip()
        meta = {}
        while not line.startswith("end type"):
          if not line == '!' and not line.startswith('! '):
            line = line.replace(", ALLOCATABLE::","")
            line = line.replace("::","")
            type, vars = re.split('\s+',line,1)
            concat = 0
            if line.startswith("type (") or line.startswith("type  ("):
              type, vars = re.split('\s+',vars,1)  # ((type.replace("type (","type(")).replace("type  (","type(")).lstrip('(')
              type = re.sub('\)','',re.sub('\(','',type))
            elif line.startswith("type("):
              type = re.sub('\)','',re.sub('type\(','',type))
            for var in re.split('[\s,]+',vars): # '\s*,(?!\:|\d)\s*',vars):
              var = var.rstrip().lstrip() # discard white space
              if var == "":
                continue
              if var.count('(') > var.count(')'):
                concat = 1
                varold = var
                continue
              if concat == 1:
                varold = varold+','+var
                if var.count('(') >= var.count(')'):
                  continue
                else:
                  concat = 0
                  var = varold
#              var = re.sub('\([0-9a-zA-Z]*,[0-9a-zA-Z]*,[0-9a-zA-Z]*\)','(:,:,:)',re.sub('\([0-9a-zA-Z]*,[0-9a-zA-Z]*\)','(:,:)',re.sub('\([0-9a-zA-Z]*\)','(:)',var)))
              if var.find('(') > 0:
                key, argum = var.split('(',1)
              else:
                key = var
              if type in listofcustomtypes.keys():
                for vart in customtypes[type].keys():
                  a, b = vart.split("~")
                  new = a+"~"+var+"~"+b
                  key2 = a+"~"+key+"~"+b
                  meta[key2] = new
              else:
                if descriptor != "":
                  tescriptions[key] = descriptor
                new = type+"~"+var
                key = type+"~"+key
                meta[key] = new
            customtypes[metatype] = meta
          else:
            line = ""
          descriptor = ""
          line = fileIN.readline()
          if line.count(" !") > 0:
            line, descriptor = line.lstrip().rstrip().split(" !",1)
            descriptor = descriptor.lstrip().rstrip()
          else:
            line = line.lstrip().rstrip()
        # move to next line
#      line = fileIN.readline().lstrip().rstrip().split(" !",10)[0]
      descriptor = ""
      line = fileIN.readline()
      if line.count(" !") > 0:
        line, descriptor = line.lstrip().rstrip().split(" !",1)
        descriptor = descriptor.lstrip().rstrip()
      else:
        line = line.lstrip().rstrip()


# loop over every mod_* file in current directory
for file in os.listdir("."):
  if (file.startswith('mod_') and not file.endswith('_interfaces.f90')):
    # open each mod_* file
    fileIN = open(file)
    line = fileIN.readline().lstrip(' ').rstrip()
    # Types (integer,etc..)
    module = {}
    # list of types temporary
    flag1 = False # add only variables inside module block
    # loop over all lines in module block
    while line:
      # build hash table of global descriptions
      if re.match("!.*?\:",line):
        type, descpt = line.split(':',1)
        type = type + " "
        x, name, y = re.split('\s+',type,2)
        descriptions[name] = descpt
      # reach end of module block
      if line.startswith("end module"):
        break
      # start adding variables
      elif line.startswith("module "):
        flag1 = True
      # cycle through type definitions
      elif line.startswith("type ") and not line.startswith("type (") and not line.startswith("type  ("):
        line = fileIN.readline().lstrip().rstrip().split(" !",10)[0]
        while line:
          if line.startswith("end type"):
            break 
          line = fileIN.readline().lstrip().rstrip()#.split(" !",10)[0]
      # remove temporary Type and change variable names to <type>%<variable>
      elif line.startswith("type(") or line.startswith("type (") or line.startswith("type  ("):
        line = (line.replace("type (","type(")).replace("type  (","type(")
        line = re.sub(",\W*ALLOCATABLE\W*::","",line)
        half, newl = re.split("\s+",line,1)
        x, metatype = half.split("(")
        metatype, x = metatype.split(")")
        dalist = customtypes[metatype].keys()
        for new in re.split("[\s,]+",newl):
          new = new+'('
          new, dummy = new.split("(",1)
          for var in dalist:
            if var.count("~") > 1:
              a, b, c = var.split("~")
              if c in tescriptions and not new+"%"+b+"%"+c in descriptions:
                descriptions[new+"%"+b+"%"+c] = tescriptions[c]
              b = new+"%"+b+"%"+c
            else:
              a, b = var.split("~")
              if b in tescriptions and not new+"%"+b in descriptions:
                descriptions[new+"%"+b] = tescriptions[b]
              b = new+"%"+b
            if module.has_key(a):
              mod = module[a]
              mod[b] = b
              module[a] = mod
              categories = {}
              can_access = {}
              access = {}
              modify = {}
              passes = {}
              passes2 = {}
              categories['can_access'] = can_access
              categories['access'] = access
              categories['modify'] = modify
              categories['passes'] = passes
              categories['passes2'] = passes2
              modulname,x = file.split('.',1)
              categories['module'] = modulname
              globals[b] = categories
            else:
              mod = {}
              mod[b] = b
              module[a] = mod
              categories = {}
              can_access = {}
              access = {}
              modify = {}
              passes = {}
              passes2 = {}
              categories['can_access'] = can_access
              categories['access'] = access
              categories['modify'] = modify
              categories['passes'] = passes
              categories['passes2'] = passes2
              modulname,x = file.split('.',1)
              categories['module'] = modulname
              globals[b] = categories
          if module.has_key(metatype):
            mod = module[metatype]
            mod[new] = new
            module[metatype] = mod
            categories = {}
            can_access = {}
            access = {}
            modify = {}
            passes = {}
            passes2 = {}
            categories['can_access'] = can_access
            categories['access'] = access
            categories['modify'] = modify
            categories['passes'] = passes
            categories['passes2'] = passes2
            modulname,x = file.split('.',1)
            categories['module'] = modulname
            globals[new] = categories
          else:
            mod = {}
            mod[new] = new
            module[metatype] = mod
            categories = {}
            can_access = {}
            access = {}
            modify = {}
            passes = {}
            passes2 = {}
            categories['can_access'] = can_access
            categories['access'] = access
            categories['modify'] = modify
            categories['passes'] = passes
            categories['passes2'] = passes2
            modulname,x = file.split('.',1)
            categories['module'] = modulname
            globals[new] = categories


      # check for parameter declarations, rewrite variables that have constant value
      elif line.startswith("parameter ("):
        x, const = line.split("(")
        const, x = const.split(")")
        param, num = const.split("=")
        dalist = module.keys()
        for var in dalist:
          mod = module[var]
          if mod.has_key(param):
            del mod[param]
            string = param+"="+num
            mod[param] = string
            module[var] = mod
            categories = {}
            can_access = {}
            access = {}
            modify = {}
            passes = {}
            passes2 = {}
            categories['can_access'] = can_access
            categories['access'] = access
            categories['modify'] = modify
            categories['passes'] = passes
            categories['passes2'] = passes2
            modulname,x = file.split('.',1)
            categories['module'] = modulname
            globals[param] = categories
            globals[string] = categories
            break
      # add all variables within module block not in a type block
      elif not line == '!' and not line.startswith('! ') and not line.startswith('#') and not line.startswith('implicit') and flag1:
        line = re.sub(",\W*ALLOCATABLE\W*::","",line)
        line = line.replace("::","")
        if not line.startswith('&'): 
          type, vars = re.split('\s+',line,1)
        else:
          line.lstrip('&')
          vars = line
        concat = 0
        for var in re.split('[\s,]+',vars): # '\s*,(?!\:|\d)\s*',vars):
          if var.count('(') > var.count(')'):
            concat = 1
            varold = var
            continue
          if concat == 1:
            varold = varold+','+var
            if var.count('(') >= var.count(')'):
              continue
            else:
              concat = 0
              var = varold
          var2 = var+"(1"
          key, x = var2.split('(',1) 
          if not type == 'data' and not type == 'DATA' and not type == 'save' and not type == 'SAVE':
            if module.has_key(type):
              mod = module[type]
              mod[key] = var
              module[type] = mod
              categories = {}
              can_access = {}
              access = {}
              modify = {}
              passes = {}
              passes2 = {}
              categories['can_access'] = can_access
              categories['access'] = access
              categories['modify'] = modify
              categories['passes'] = passes
              categories['passes2'] = passes2
              modulname,x = file.split('.',1)
              categories['module'] = modulname
              globals[key] = categories
              globals[var] = categories
            else:
              mod = {}
              mod[key] = var
              module[type] = mod
              categories = {}
              can_access = {}
              access = {}
              modify = {}
              passes = {}
              passes2 = {}
              categories['can_access'] = can_access
              categories['access'] = access
              categories['modify'] = modify
              categories['passes'] = passes
              categories['passes2'] = passes2
              modulname,x = file.split('.',1)
              categories['module'] = modulname
              globals[key] = categories
              globals[var] = categories
      # move to next line
      descriptor = ""
      line = fileIN.readline()
      if line.count(" !") > 0:
        line, descriptor = line.lstrip().rstrip().split(" !",1)
        descriptor = descriptor.lstrip().rstrip()
      else:
        line = line.lstrip().rstrip()

    # remove temporary keys    
    for type in listoftypes.keys():
      del module[type]   
    # combine tables for each file
    modules[file] = module

#################################################################################
#
#  Extension to script that goes through all routines (*.f90) and defines what
#  globals that routine 'can access', 'does access', and 'modifies.'
#
#################################################################################



# can access
files = {}
routine = {}
for file in os.listdir("."):
  if file.endswith('.f90') and not file.startswith('mod_'):

    subroutines = {}
    fileIN = open(file)
    notnewflag = False

    for line in fileIN:
      line = line.lstrip(' ').rstrip()
      # enter old can_access and select new routine name
      if line.startswith("subroutine ") or line.startswith("program ") or line.startswith("function "):
        if notnewflag:
          routine[routinename] = can_access
        notnewflag = True
        xtp, routinename = re.split('\s+',line,1)
        routinename = routinename + '('
        routinename, x = routinename.split('(',1)
        subroutines[routinename] = xtp+'_'+routinename
        can_access = {}
      elif line.startswith("use "):
        a,b = re.split(' | !|\!',line)
        b = "mod_" + b + ".f90"
        if modules.has_key(b):
          mod = modules[b]
          dalist = mod.keys()
          for key in dalist:
            if not can_access.has_key(key):
              typetable = {}
              can_access[key] = typetable
            module = mod[key]
            varlist = module.keys()
            for var in varlist:
              value = module[var]
              typetable = can_access[key]
              var = var.replace('(','').replace(')','').replace(':','').replace(',','')
              var = var.lstrip().rstrip()
              if var == "":
                continue
              typetable[var] = value
              can_access[key] = typetable
    routine[routinename] = can_access
    files[file] = subroutines
#print routine
#print ">>>>>>>>>>>>>>>>>>>>>"
# Print table ######################
'''routinelist = routine.keys()
for r in routinelist:
  print r
  typelist = routine[r].keys()
  for t in typelist:
    print "   <"+t+">"
    types = routine[r]
    vars = ""
    varlist = types[t].keys()
    for v in varlist:
      vars = vars + " " + v
    print "     "+vars+"\n"
'''# make list of accessed
#####################################
#print "<<<<<<<<<<<<<<<<<<<<< can access"

routine_access = {}
for file in os.listdir("."):
  if file.endswith('.f90') and not file.startswith('mod_'):

    fileIN = open(file)
    notnewflag = False

    for line in fileIN:
      line = line.lstrip(' ').rstrip()
      if line.startswith("subroutine ") or line.startswith("program ") or line.startswith("function "):
        if notnewflag:
          routine_access[routinename] = access
        notnewflag = True
        x, routinename = re.split('\s+',line,1)
        routinename = routinename + '('
        routinename, x = routinename.split('(',1)
        access = {}
        types = routine[routinename]
      elif not (line=='!') and not line.startswith('! ') and notnewflag:
        list = re.split('\W+',line)
        listcomb = []
        pos3 = ''
        for p1 in range(1,len(list)-1):
          pos1 = list[p1]
          for p2 in range(p1+1,len(list)):
            pos2 = list[p2]
            pos3 = pos1 + '%' +pos2
            listcomb.append(pos3)
            for p4 in range(p2+1,len(list)):
              pos4 = list[p4]
              pos5 = pos1 + '%' +pos2 + '%' +pos4
              listcomb.append(pos5)
        list.extend(listcomb)
        typelist = types.keys()
        for each in list:
          each = each.lstrip().rstrip()
          if each == "":
            continue
          for type in typelist:
            vars = types[type]
            if vars.has_key(each):
              if not access.has_key(type):
                access[type] = {}
              value = vars[each]
              temp1 = globals[value]
              temp2 = temp1['access']
              temp2[routinename] = routinename
              temp1['access'] = temp2
              globals[value] = temp1
              access_type = access[type]
              access_type[each] = value
              access[type] = access_type
    routine_access[routinename] = access
#######################################
'''routinelist = routine_access.keys()
for r in routinelist:
  print r
  typelist = routine_access[r].keys()
  for t in typelist:
    print "   <"+t+">"
    types = routine_access[r]
    vars = ""
    varlist = types[t].keys()
    for v in varlist:
      vars = vars + " " + v
    print "     "+vars+"\n"
'''#######################################


routine_passes = {}
routine_calls = {}
for file in os.listdir("."):
  if file.endswith('.f90') and not file.startswith('mod_'):

    fileIN = open(file)
    notnewflag = False

    for line in fileIN:
      line = line.lstrip(' ').rstrip()
      if line.startswith("subroutine ") or line.startswith("program ") or line.startswith("function "):
        if notnewflag:
          routine_passes[routinename] = passes
          routine_calls[routinename] = calls
        notnewflag = True
        x, routinename = re.split('\s+',line,1)
        routinename = routinename + '('
        routinename, x = routinename.split('(',1)
        passes = {}
        calls = []
        types = routine[routinename]
      elif not (line=='!') and not line.startswith('! ') and notnewflag:
        list = re.split('\W+',line)
        listcomb = []
        if ("call" in list) or ("Call" in list) or ("CALL" in list):
          callstr = ""
          if ("call" in list):
            callstr = "call"
          if ("Call" in list):
            callstr = "Call"
          if ("CALL" in list):
            callstr = "CALL"
          for p1 in range(list.index(callstr)+1,len(list)-1):
            pos1 = list[p1]
            for p2 in range(p1+1,len(list)):
              pos2 = list[p2]
              pos3 = pos1 + '%' +pos2
              listcomb.append(pos3)
              for p4 in range(p2+1,len(list)):
                pos4 = list[p4]
                pos5 = pos1 + '%' +pos2 + '%' +pos4
                listcomb.append(pos5)
            
          list.extend(listcomb)
          typelist = types.keys()
          for p1 in range(list.index(callstr)+1,len(list)):
            each = list[p1]
            each = each.lstrip().rstrip()
            if each == "":
              continue
            for type in typelist:
              vars = types[type]
              if vars.has_key(each):
                if not passes.has_key(type):
                  passes[type] = {}
                
                value = vars[each]
                temp1 = globals[value]
                temp2 = temp1['passes']
                temp2[routinename] = routinename
                temp1['passes'] = temp2
                temp2 = temp1['passes2']
                temp2[list[list.index(callstr)+1]] = list[list.index(callstr)+1]
                temp1['passes2'] = temp2
                globals[value] = temp1
                 
                passes_type = passes[type]
                passes_type[each] = value
                passes[type] = passes_type
                if not list[list.index(callstr)+1] in calls:
                  calls.append(list[list.index(callstr)+1])

    routine_passes[routinename] = passes
    routine_calls[routinename] = calls



routine_modify = {}
for file in os.listdir("."):
  if file.endswith('.f90') and not file.startswith('mod_'):

    fileIN = open(file)
    notnewflag = False

    for line in fileIN:
      line = line.lstrip(' ').rstrip()
      if line.startswith("subroutine ") or line.startswith("program ") or line.startswith("function "):
        if notnewflag:
          routine_modify[routinename] = modify
        notnewflag = True
        x, routinename = re.split('\s+',line,1)
        routinename = routinename + '('
        routinename, x = routinename.split('(',1)
        modify = {}
        types = routine_access[routinename]
      elif not (line=='!') and not line.startswith('! ') and (line.find('=')>=0):
        list = re.split('=',line)
        list.pop()
        for each in list:
          each = each.rstrip()
          if each.endswith("'"):
            continue
          each = re.sub('\([a-zA-Z0-9,:%+-]*?\)','',each)
          each = re.sub('\([a-zA-Z0-9,:%+-]*?\)','',each)
          each = re.sub('\([a-zA-Z0-9,:%+-]*?\)','',each)
          each = re.sub('\([a-zA-Z0-9,:%+-]*?\)','',each)
#          listi = each + '('
#          listi, x = listi.split('(',1)
#          listi = re.split('\W+',listi)
          list2 = re.split('[^a-zA-Z0-9%_]+',each)
          listcomb = []
          pos3 = ''
          for p1 in range(1,1):
            pos1 = list2[p1]
            for p2 in range(p1+1,len(list2)):
              pos2 = list2[p2]
              pos3 = pos1 + '%' +pos2
              listcomb.append(pos3)
              for p4 in range(p2+1,len(list2)):
                pos4 = list2[p4]
                pos5 = pos1 + '%' +pos2 + '%' +pos4
                listcomb.append(pos5)
#          listi.extend(listcomb)
#          list2 = listi
          typelist = types.keys()
          for each2 in list2:
            each2 = each2.lstrip().rstrip()
            if each2 == "":
              continue
            for type in typelist:
              vars = types[type]
              if vars.has_key(each2):
                if not modify.has_key(type):
                  modify[type] = {}
                value = vars[each2]
                temp1 = globals[value]
                temp2 = temp1['modify']
                temp2[routinename] = routinename
                temp1['modify'] = temp2
                globals[value] = temp1
                modify_type = modify[type]
                modify_type[each2] = value
                modify[type] = modify_type
    routine_modify[routinename] = modify
#print ">>>>>>>>>>>>>>>>>>>>>>>>>"
#print routine_modify

# Print globals and where they are accessed or modified ####
'''globs = globals.keys()
for g in globs:
  print g
  print "  access"
  temp1 = globals[g]
  types = temp1['access'].keys()
  for t in types:
    print "    " + t
  print "  modify"
  types = temp1['modify'].keys()
  for t in types:
    print "    " + t
'''#
###########################################################

# Print table ######################
'''filelist = files.keys()
for file in filelist:
  print file
  routinelist = files[file].keys()
  for r in routinelist:
    print "  ________ " +r+" ________"
    print "    {access}"
    typelist = routine_access[r].keys()
    for t in typelist:
      types = routine_access[r]
      vars = "        <"+t+">"
      varlist = types[t].keys()
      for v in varlist:
        vars = vars + "; " + v
      print vars
    print "    {modify}"
    typelist = routine_modify[r].keys()
    for t in typelist:
      vars =  "       <"+t+">"
      types = routine_modify[r]
      varlist = types[t].keys()
      for v in varlist:
        vars = vars + "; " + v
      print vars
  print "############################################"
'''# make list of accessed
#####################################
#print module files

for i in globals.keys():
  if (i.find('%')>=0):
    i45, i44 = i.split('%',1)
    i45 = i45 + "(1"
    i43, x = i45.split('(',1)
    i44 = i44 + "(1"
    i42, x = i44.split('(',1)
    i3 = i43 + '%' + i42 +"(1"
  else:
    i3 = i + "(1"
  i2, x = i3.split('(',1)
  globals[i2] = globals[i]

moduleslist = modules.keys()
for module in moduleslist:
  if not os.path.isdir('globals_doc'):
    os.mkdir('globals_doc')
  modulename, x = module.split('.',1)
  filehandle = open('globals_doc/'+modulename+'.html','w')
  filehandle.write('<html>')
  filehandle.write('<head></head>')
  filehandle.write('<body>')
  filehandle.write('<table width="100%" border="1" cellpadding="3" cellspacing="0">')
  filehandle.write('<tbody><tr class TableHeadingColor" bgcolor="#cccfff">')
  filehandle.write('<th colspan="1" align="left"><font size="+2"><b>'+module)
  filehandle.write('</b></font></th></tr></tbody></table>')
  filehandle.write('<a href="global_tracker.html">(index)</a>')
  mod = modules[module]
  typelist = mod.keys()
  for t in typelist:
    if not t == 'data' and not type == 'DATA' and not type == 'save' and not type == 'SAVE':
      typelist = mod[t]
      vars = ""
      varlist = typelist.keys()
      varlist.sort()
#      filehandle.write('<dl>')
      for v in varlist:
        vars = vars + "; " + typelist[v]
        dtt = typelist[v]
        dtt = dtt + "=1"
        typestrip, x = dtt.split('=',1)
        if not (typestrip.find('%')>=0):
          typestrip = typestrip + '('
          typestrip, x = typestrip.split('(',1)
        typestrip = typestrip.replace('(','').replace(')','').replace(':','').replace(',','')
        if  descriptions.has_key(typestrip):
          descptkey = descriptions[typestrip]
        else:
          descptkey = ' '
        filehandle.write('<dt><a name="'+typestrip+'"></a><h3>' + '<span style="color:green">'+t+'  </span>'+typelist[v] +'</h3></dt><dd><dl>')
        filehandle.write('<dt>'+descptkey+'</dt>')
        filehandle.write('<dt><b>Accessed by:</b></dt><dd>')
        methodlist = globals[typelist[v]]
        sublist = methodlist['access'].keys()
        sublist.sort()
        for s in sublist:
          filename = modulename
          filelist = files.keys()
          for file in filelist:
            if files[file].has_key(s):
              filename, x = file.split('.',1)
          filehandle.write('<a href="'+filename+'.html#'+s+'">'+s+'</a>, ')
        filehandle.write('</dd>')
        filehandle.write('<dt><b>Passed by:</b></dt><dd>')
        methodlist = globals[typelist[v]]
        sublist = methodlist['passes'].keys()
        sublist.sort()
        for s in sublist:
          filename = modulename
          filelist = files.keys()
          for file in filelist:
            if files[file].has_key(s):
              filename, x = file.split('.',1)
          filehandle.write('<a href="'+filename+'.html#'+s+'">'+s+'</a>, ')
        filehandle.write('</dd>')
        filehandle.write('<dt><b>Passed to:</b></dt><dd>')
        methodlist = globals[typelist[v]]
        sublist = methodlist['passes2'].keys()
        sublist.sort()
        for s in sublist:
          filename = ""
          filelist = files.keys()
          for file in filelist:
            if files[file].has_key(s):
              filename, x = file.split('.',1)
          if (filename == ""):
            filehandle.write(s+' ')
          else:
            filehandle.write('<a href="'+filename+'.html#'+s+'">'+s+'</a>, ')
        filehandle.write('</dd>')
        filehandle.write('<dt><b>Directly modified by:</b></dt><dd>')
        methodlist = globals[typelist[v]]
        sublist = methodlist['modify'].keys()
        sublist.sort()
        for s in sublist:
          filename = modulename
          filelist = files.keys()
          for file in filelist:
            if files[file].has_key(s):
              filename, x = file.split('.',1)
          filehandle.write('<a href="'+filename+'.html#'+s+'">'+s+'</a>, ')
        filehandle.write('</dd>')
        filehandle.write('</dd></dd></dl><dt><hr></dt>')
#      filehandle.write('</dl>')
  filehandle.write('</body>')
  filehandle.write('</html>')

#####################################
#print routine files
filelist = files.keys()
for file in filelist:
  if not os.path.isdir('globals_doc'):
    os.mkdir('globals_doc')
  filename, x = file.split('.',1)
  filehandle = open('globals_doc/'+filename+'.html','w')
  filehandle.write('<html>')
  filehandle.write('<head></head>')
  filehandle.write('<body>')
  filehandle.write('<table width="100%" border="1" cellpadding="3" cellspacing="0">')
  filehandle.write('<tbody><tr class TableHeadingColor" bgcolor="#fffccc">')
  filehandle.write('<th colspan="1" align="left"><font size="+2"><b>'+file)
  filehandle.write('</b></font></th></tr></tbody></table>')
  filehandle.write('<a href="global_tracker.html">(index)</a>')
  routinelist = files[file].values()
  for rr in routinelist:
    subword,r = re.split('_',rr,1)
    filehandle.write('<dl><dt><a name="'+r+'"></a><h3>' + '<span style="color:purple">'+subword+'  </span>'+ r +'</h3></dt><dd><dl>')
    filehandle.write('<dt><b>Access:</b></dt>')
    typelist = routine_access[r].keys()
    for t in typelist:
      filehandle.write('<dd>{'+t+'} ')
      typeg = routine_access[r]
      varlist = typeg[t].keys()
      for v in varlist:
        gb = globals[v]
        filehandle.write('<a href="'+gb['module']+'.html#'+v+'">'+v+'</a>, ')
    filehandle.write('<dt><b>Modify:</b></dt>')
    typelist = routine_modify[r].keys()
    for t in typelist:
      filehandle.write('<dd>{'+t+'} ')
      typeg = routine_modify[r]
      varlist = typeg[t].keys()
      for v in varlist:
        gb = globals[v]
        filehandle.write('<a href="'+gb['module']+'.html#'+v+'">'+v+'</a>, ')
    filehandle.write('<dt><b>Passes:</b></dt>')
    typelist = routine_passes[r].keys()
    for t in typelist:
      typeg = routine_passes[r]
      varlist = typeg[t].keys()
      for v in varlist:
        gb = globals[v]
        p2l = (gb['passes2']).keys()
        filelist = files.keys()
        filename, xx = file.split('.',1)
#        for file2 in filelist:
#          if files[file2].has_key(''.join(p2l)):
#            filename, xx = file2.split('.',1)
        filehandle.write('<dd>{'+t+'} <a href="'+gb['module']+'.html#'+v+'">'+v+'</a> to ')
        for subr in p2l:
          if not subr in routine_calls[r]:
            continue
          filename2 = ""
          for file2 in filelist:
            if subr in files[file2].keys():
              filename2, xx = file2.split('.',1)
          if filename2 == "":
            filehandle.write(subr+', ')
          else:
            filehandle.write('<a href="'+filename2+'.html#'+subr+'">'+' '+subr+'</a>, ')
    filehandle.write('</dd></dl></dd></dl><hr>')
  filehandle.write('</body>')
  filehandle.write('</html>')

#####################################
# make index file
if not os.path.isdir('globals_doc'):
  os.mkdir('globals_doc')
filehandle = open('globals_doc/global_tracker.html','w')
filehandle.write('<html>')
filehandle.write('<head></head>')
filehandle.write('<body>')
filehandle.write('<h1> CAMPARI Code - Global Variable Tracker </h1>')
filehandle.write('<dl><dt><h3> Subroutine File Index: </h3></dt>')
filelist = files.keys()
filelist.sort()
for file in filelist:
  filename, x = file.split('.',1)
  filehandle.write('<dd><a href="'+filename+'.html">'+file+'</a></dd>')
filehandle.write('<dt><h3> Module File Index: </h3></dt>')
modulelist = modules.keys()
modulelist.sort()
for module in modulelist:
  modulename, x = module.split('.',1)
  filehandle.write('<dd><a href="'+modulename+'.html">'+module+'</a></dd>')
filehandle.write('</dl>')
filehandle.write('</body>')
filehandle.write('</html>')
