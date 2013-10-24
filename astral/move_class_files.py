import sys, codecs, os, re
if sys.stdout.encoding is None:
        sys.stdout = codecs.open('/dev/stdout', 'w', 'utf-8')

profiles_dir = '../spherical_system/profiles/'
profiles = os.listdir(profiles_dir)

fr = open('class-c','r')
num_lines = 0
count = 0
line = fr.readline()
while line != '':
  x = line.strip('\n')
  y = line.split()
  domain = y[0]
  file_name = y[0] + '.ent.profile'
  if file_name in profiles:
    count += 1
    cmd = 'cp ' + profiles_dir + file_name + ' ../spherical_system/class-c/'
    os.system(cmd)
  line = fr.readline()
  num_lines += 1

print '# of lines: ', num_lines
print 'count: ', count    
