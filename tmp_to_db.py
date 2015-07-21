columns = ["color", "chr", "v", "k", "lambda", "mu", "r_power_f", "s_power_g", "comments"]

with open('brouwer.tmp') as f:
    lines = f.readlines()

data = [dict(zip(columns,l.strip().split('|'))) for l in lines]

color_to_status = {'#b0ffb0':'exists', '#ffb0b0':'impossible', '#ffffb0':'open'}

for i,x in enumerate(data):
    x['status'] = color_to_status[x.pop('color')]
    if x['v'] == '':
        x['v'] = data[i-1]['v']
    for key in ['v','k','lambda','mu']:
        x[key] = int(x[key])

data = {(x.pop('v'),x.pop('k'),x.pop('lambda'),x.pop('mu')):x
        for x in data}

# Text output
with open('brouwer.txt','w') as output:
    for (v,k,l,mu),dic in sorted(data.items()):
        output.write("{:<4} {:<4} {:<4} {:<4} {:<10} {}\n".format(v,k,l,mu,dic['status'],dic['comments']))

print "'brouwer.txt' file written."

# sobj output
from sage.structure.sage_object import save
save(data,'brouwer.sobj')

print "'brouwer.sobj' file written."

stats = [x['status'] for x in data.values()]
print "statistics:"
print " - {} impossible".format(stats.count('impossible'))
print " - {} open".format(stats.count('open'))
print " - {} realizable".format(stats.count('exists'))
