#!/usr/bin/python
__author__="alvaro"
__date__ ="$Sep 2, 2010 10:09:19 AM$"

import sys;
import re;
import shlex;
import os;
from os import path;
from optparse import OptionParser;
from subs import subs;

KEYWORD = '@precisions';
REGEX = '^.*'+KEYWORD+'\s+((\w+,?)+)\s+(\w+)\s+->\s*((\s\w+)+).*$';

def relpath(p):
  p = path.realpath(p);
  return p.replace(path.realpath('.')+'/','');

class Conversion:
  debug = False;
  test = False;
  make = False;
  prefix = None;
  required_precisions = [];
  files_in = [];
  files_out = [];
  def __init__(self, file, match, content):
    self.content = content;
    file = path.realpath(file);
    rel = relpath(file);
    self.file = list(path.split(file));
    self.date = path.getmtime(file);
    if path.samefile(path.join(self.file[0],self.file[1]),sys.argv[0]):
      raise ValueError('Let\'s just forget codegen.py');
    try:
      # normal,all,mixed
      self.types = match[0].split(',');
      # z
      self.precision = match[2].lower();
      # ['c','d','s']
      self.precisions = match[3].lower().split();
      if len(self.required_precisions):
        self.precstmp = [];
        for prec in (self.required_precisions):
          if prec in(self.precisions):
            self.precstmp.append(prec);
        self.precisions = self.precstmp;
    except:
      raise ValueError(path.join(self.file[0],self.file[1])+' : Invalid conversion string');
    self.files_in.append(rel);


  def run(self):
    if self.convert_names() and not self.test:
      self.convert_data();
      self.export_data();

  def convert_names(self):
    self.names = [];
    self.dates = [];
    self.converted = [];
    load = False;
    if self.debug: print '|'.join(self.types), self.precision, relpath(path.join(self.file[0],self.file[1]));
    for precision in self.precisions:
      new_file = self.convert(self.file[1], precision);
      if self.debug: print precision,':',
      if new_file <> self.file[1]:
        if self.prefix is None:
          prefix = self.file[0]
          makeprefix = '';
        else:
          prefix = self.prefix;
          makeprefix = '--prefix '+prefix;
        conversion = path.join(prefix, new_file);
        file_out = relpath(conversion);
        if self.make:
          file_in = relpath(path.join(self.file[0],self.file[1]));
          print file_out+':',file_in;
          print "\t$(PYTHON)",path.realpath(sys.argv[0]),makeprefix,'-p',precision,"--file",file_in;
        self.names.append(new_file);
        self.files_out.append(file_out);
        if self.debug: print relpath(conversion), ':',
        try:
          date = path.getmtime(conversion);
          diff = self.date - date;
          self.dates.append(diff);
          if self.debug:
            if diff > 0: print 'Old',
            else: print 'Current',
            print diff;
          if diff > 0: load = True;
        except:
          if self.debug: print 'Missing';
          self.dates.append(None);
          load = True;
      else:
        if self.debug: print '<No Change>',':';
        else: print >> sys.stderr, new_file, 'had no change for', precision;
        self.names.append(None);
        self.dates.append(None);
    return load;

  def export_data(self):
    for i in range(len(self.names)):
      name = self.names[i];
      data = self.converted[i];
      if data is None or name is None: continue;
      if self.prefix is None:
        fd = open(path.join(self.file[0],name), 'w');
      else:
        fd = open(path.join(self.prefix,name), 'w');
      fd.write(data);
      fd.close();


  def convert_data(self):
    for i in range(len(self.precisions)):
      precision = self.precisions[i];
      name = self.names[i];
      date = self.dates[i];
      if name is not None and (date is None or date > 0):
        self.converted.append(self.convert(self.content, precision));
      else: self.converted.append(None);

  def substitute(self, sub_type, data, precision):
    try:
      work = subs[sub_type];
      prec_to = work[0].index(precision);
      prec_from = work[0].index(self.precision);
    except:
      return data;
    for i in range(1,len(work)):
      try:
        search = work[i][prec_from];
        replace = work[i][prec_to];
        if not search: continue;
        replace.replace('\*','*');
        data = re.sub(search, replace, data);
      except:
        print 'Bad replacement pair ',i,'in',sub_type;
        continue;
    return data;

  def convert(self, data, precision):
    try:
      data = self.substitute('all', data, precision);
    except: pass;
    for sub_type in self.types:
      if sub_type == 'all': continue;
      try:
        data = self.substitute(sub_type, data, precision);
      except Exception, e:
        raise ValueError('I encountered an unrecoverable error while working in subtype:',sub_type+'.');
    data = re.sub('@precisions '+','.join(self.types)+'.*', '@generated '+precision, data); 
    return data;

def grep(string,list):
    expr = re.compile(string)
    return filter(expr.search,list)

parser = OptionParser();
parser.add_option("-d", "--debug",     help='Print debugging messages.',                                   action='store_true', dest='debug',     default=False);
parser.add_option("-P", "--prefix",    help='The output directory if different from the input directory.', action='store',      dest='prefix',    default=None);
parser.add_option("-i", "--in-files",  help='Print the filenames of files for precision generation.',      action='store_true', dest='in_print',  default=False);
parser.add_option("-o", "--out-files", help='Print the filenames for the precision generated files.',      action='store_true', dest='out_print', default=False);
parser.add_option("-c", "--clean",     help='Remove the files that are the product of generation.',        action='store_true', dest='out_clean', default=False);
parser.add_option("-t", "--threads",   help='Enter the number of threads to use for conversion.',          action='store',      dest='threads',   default=1);
parser.add_option("-f", "--file",      help='Specify a file(s) on which to operate.',                      action='store',      dest='fileslst', type='string', default="");
parser.add_option("-p", "--prec",      help='Specify a precision(s) on which to operate.',                 action='store',      dest='precslst', type='string', default="");
parser.add_option("-m", "--make",      help='Spew a GNU Make friendly file to standard out.',              action='store_true', dest='make',      default=False);
parser.add_option("-T", "--test",      help='Don\'t actually do any work.',                                action='store_true', dest='test',      default=False);
(options, args) = parser.parse_args();

rex = re.compile(REGEX);
work = [];

def check_gen(file):
  fd = open(path.realpath(file), 'r');
  lines = fd.readlines();
  fd.close();
  for line in lines:
    m = rex.match(line);
    if m is None: continue;
    work.append((file, m.groups(), ''.join(lines)));

if options.fileslst:
  options.files = options.fileslst.split();
  if len(options.files):
    for file in options.files:
      check_gen(file);
else:
  startDir = '.';
  for root, dirs, files in os.walk(startDir, True, None):
    for file in files:
      if file.startswith('.'): continue;
      if not file.endswith('.c') and not file.endswith('.h') and not file.endswith('.f'):
        continue;
      check_gen(path.join(root,file));
    if '.svn' in dirs:
      dirs.remove('.svn');

Conversion.debug = options.debug;
Conversion.make = options.make;
Conversion.prefix = options.prefix;

if options.out_print or options.out_clean or options.in_print or options.make or options.test:
  Conversion.test = True;

if options.precslst:
  options.precs = options.precslst.split();
  if len(options.precs):
    Conversion.required_precisions = options.precs;
if options.make:
  print '## Automatically generated Makefile';
  print 'PYTHON ?= python';

if len(work) is 0:
  print 'gen = ';
  print 'cleangen:';
  print '\trm -f $(gen)';
  print 'generate: $(gen)';
  print '.PHONY: cleangen generate';
  sys.exit(0);

for tuple in work:
  try:
    c = Conversion(tuple[0], tuple[1], tuple[2]);
    c.run();
  except Exception, e:
    print >> sys.stderr, str(e);
    continue;

if options.make:
  print 'gen = ',' '+' '.join(c.files_out);
  print 'cleangen:';
  print '\trm -f $(gen)';
  print 'generate: $(gen)';
  print '.PHONY: cleangen generate';
if options.in_print: print ' '.join(c.files_in);
if options.out_print: print ' '.join(c.files_out);
if options.out_clean:
  for file in c.files_out:
    if not path.exists(file): continue;
    os.remove(file);
