import sys, codecs, os, re
if sys.stdout.encoding is None:
        sys.stdout = codecs.open('/dev/stdout', 'w', 'utf-8')

fr = open('astral-scopdom-seqres-gd-sel-gs-bib-40-1.75B.fa.txt','r')
fw = open('domains.test','w')
line = fr.readline()
num_lines = 0
num_domains = 0
classes = {}

while line != '':
  x = line.strip('\n')
  y = line.split()
  if (y[0][0] == '>'):
    num_domains += 1
    domain = y[0][1:]
    domain_id = y[1]
    class_id = y[1].split('.')[0]
    if class_id in classes:
      classes[class_id].append(domain)
    else:
      classes[class_id] = [domain]
    fw.write(domain+'\t'+domain_id+'\n')
  num_lines += 1
  line = fr.readline()

fw.close()
fr.close()
print '# of lines: ', num_lines
print '# of domains: ', num_domains
print '# of classes: ', len(classes)
for each_class in classes:
  print each_class, len(classes[each_class])
  file_name = 'class-' + each_class 
  fw = open(file_name,'w')
  for domain in classes[each_class]:
    fw.write(domain+'\n')
  fw.close()

