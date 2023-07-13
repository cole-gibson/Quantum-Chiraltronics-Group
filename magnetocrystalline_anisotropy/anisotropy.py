import numpy as np
from sklearn.linear_model import LinearRegression
import scipy
from tabulate import tabulate
# Set path to file with energies (theta, phi, energy)
ENERGYPATH = 'energies'
INPUT = 'input'

Dict = dict()

# How to handle multiple qs. Treated as one term with one coefficient

for i in ['Oh', 'O', 'Th', 'Td', 'T']:
    Dict[i] = ['a1**2*a2**2 + a2**2*a3**2 + a1**2*a3**2']
    
for i in ['Oh', 'O', 'Th']:
    Dict[i].append('a1**2*a2**2*a3**2')
    
for i in ['Td', 'T']:
    Dict[i].append('a1*a2*a3')
    
for i in ['Th', 'T']:
    Dict[i].append('a1**4*(a2**2-a3**2) + a2**4*(a3**2-a1**2) + a3**4*(a1**2 - a2**2)')
    
Dict['O'].append('a1*a2*a3*(a1**4*(a2**2-a3**2) + a2**4*(a3**2-a1**2) + a3**4*(a1**2 - a2**2))')

for i in ['D4h', 'D2d', 'C4h', 'D4', 'S4', 'D6h', 'D3d', 'C6h', 'D6', 'S6']:
    Dict[i] = ['a3**2']
    
for i in ['D4h', 'D2d', 'C4h', 'D4', 'S4']:
    Dict[i].append('a1**2*a2**2')

for i in ['D6h', 'D3d', 'C6h', 'D6', 'S6']:
    Dict[i].append('a1**6 - 15*a1**4*a2**2 + 15*a1**2*a2**4 - a2**6')
    
Dict['D2d'].append('a1*a2*a3')
Dict['C4h'].append('a1*a2*(a1**2 - a2**2)')
Dict['D4'].append('a1*a2*a3*(a1**2 - a2**2)')
Dict['S4'].append('a1*a2*a3 + a3*(a1**2 - a2**2) + a1*a2*(a1 - a2**2)') # mistake here?

Dict['D3d'].append('a3*(a1**3 - 3*a1*a2**2)')
Dict['C6h'].append('a1*a2*(3*a1**4 - 10*a1**2*a2**2 + 3*a2**4)')
Dict['D6'].append('a1*a2*a3*(3*a1**4 - 10*a1**2*a2**2 + 3*a2**4)')
Dict['S6'].append('a3*(a1**3 - 3*a1*a2**2) + a3*(a2**3 - 3*a1**2*a2) + (a1**3 - 3*a1*a2**2)*(a2**3 - 3*a1**2*a2)')

for i in ['C6v', 'C6']:
    Dict[i] = ['a3']
    Dict[i].append('a1**6 - 15*a1**4*a2**2 + 15*a1**2*a2**4 - a2**6')
    
Dict['C6'].append('a1*a2*(3*a1**4 - 10*a1**2*a2**2 + 3*a2**4')

for i in ['D2h', 'C2h', 'D2', 'Ci']:
    Dict[i] = ['a1**2']
    Dict[i].append('a2**2')
    
Dict['C2h'].append('a1*a2')
Dict['D2'].append('a1*a2*a3')
Dict['Ci'].append('a1*a2 + a2*a3 + a3*a1')

for i in ['C4v', 'C4']:
    Dict[i] = ['a3']
    Dict[i].append('a1**2*a2**2')
    
Dict['C4'].append('a1*a2*(a1**2 - a2**2)')

for i in ['D3h', 'C3h', 'D3']:
    Dict[i] = ['a3**2']
    Dict[i].append('a1**3 - 3*a1*a2**2')
    
Dict['C3h'].append('a2**3 - 3*a1**2*a2')
Dict['D3'].append('a3*(a2**3 - 3*a1**2*a2)')

for i in ['C3v', 'C3']:
    Dict[i] = ['a3']
    Dict[i].append('a1**3 - 3*a1*a2**2')

Dict['C3'].append('a2**3 - 3*a1**2*a2')

for i in ['C2v', 'C2']:
    Dict[i] = ['a3']
    Dict[i].append('a1**2 - a2**2')
    
Dict['C2'].append('a1*a2')

for i in ['Cs', 'C1']:
    Dict[i] = ['a1']
    Dict[i].append('a2')

Dict['C1'].append('a3')

with open(INPUT, 'r') as f:
    lines = f.readlines()
    symmetry = lines[0].replace('\n', '')
    for i in range(1, len(lines)):
        if lines[i].strip():
            break
    d1 = eval(lines[i].replace('\n', ''))
    d2 = eval(lines[i+1].replace('\n', ''))
    d3 = eval(lines[i+2].replace('\n', ''))

energy = np.loadtxt(ENERGYPATH)

Doring = Dict

termsin = Doring[symmetry].copy()

for i in range(len(termsin)):
    for j in range(3):
        termsin[i] = termsin[i].replace('a{}'.format(j+1), 'cos[:,{}]'.format(j))

tolterms = len(termsin)

def normalize(x):
    return np.divide(x,np.linalg.norm(x))

axis = np.array([d1,d2,d3])

for i in range(len(axis)):
    axis[i] = normalize(axis[i])

l,w=energy.shape

positions=np.zeros((l,3))
positions[:,0]=np.cos(energy[:,1])*np.sin(energy[:,0])
positions[:,1]=np.sin(energy[:,1])*np.sin(energy[:,0])
positions[:,2]=np.cos(energy[:,0])

cos=np.zeros((l,3))
for i in np.arange(0,l):
    positions[i,:]=normalize(positions[i,:])
    for j in np.arange(0,3):        
        cos[i,j]=np.dot(positions[i,:],axis[j])

terms=np.zeros((l, tolterms))
for i in range(tolterms):
    terms[:,i] = eval(termsin[i])

# The fitting
reg=LinearRegression().fit(terms,energy[:,2])
reg.score(terms,energy[:,2])

# Print fitting equation
cosplot=np.array(["nan","nan","nan"],dtype="object")
point=np.array(["cos(v)*sin(u)","sin(v)*sin(u)","cos(u)"])
for j in np.arange(0,3):        
    cosplot[j]=str(axis[j,0])+"*"+point[0]+"+"+str(axis[j,1])+"*"+point[1]+"+"+str(axis[j,2])+"*"+point[2]

def x0(u,v):
    return axis[0,0]*np.cos(v)*np.sin(u)+axis[0,1]*np.sin(v)*np.sin(u)+axis[0,2]*np.cos(u)

def x1(u,v):
    return axis[1,0]*np.cos(v)*np.sin(u)+axis[1,1]*np.sin(v)*np.sin(u)+axis[1,2]*np.cos(u)

def x2(u,v):
    return axis[2,0]*np.cos(v)*np.sin(u)+axis[2,1]*np.sin(v)*np.sin(u)+axis[2,2]*np.cos(u)

def f(params):
    u,v=params
    return reg.intercept_+reg.coef_[0]*((x0(u,v))**2*(x1(u,v))**2+(x0(u,v))**2*(x2(u,v))**2+(x1(u,v))**2*(x2(u,v))**2)+reg.coef_[1]*(x0(u,v))**2*(x1(u,v))**2*(x2(u,v))**2

def negf(params):
    return -1*f(params)

# Find min and max values of energy
minimum=scipy.optimize.minimize(f,[1,1],bounds=((0,np.pi),(0,2*np.pi)))
u,v=minimum.x
easy=[x0(u,v),x1(u,v),x2(u,v)]
easy = easy/max(list(map(abs, easy)))

maximum=scipy.optimize.minimize(negf,[1,1],bounds=((0,np.pi),(0,2*np.pi)))
u,v=maximum.x
hard=[x0(u,v),x1(u,v),x2(u,v)]
hard = hard/max(list(map(abs, hard)))

termplot = termsin

for i in range(len(termplot)):
    for j in range(3):
        termplot[i] = termplot[i].replace('cos[:,{}]'.format(j), '('+str(cosplot[j])+')')
# Old
# print(str(-1*min.fun)+"+"+str(reg.intercept_)+"+"+str(reg.coef_[0])+"*(("+cosplot[0]+")**2*("+cosplot[1]+")**2+("+cosplot[0]+")**2*("+cosplot[2]+")**2+("+cosplot[1]+")**2*("+cosplot[2]+")**2)"+"+"+str(reg.coef_[1])+"*("+cosplot[0]+")**2*("+cosplot[1]+")**2*("+cosplot[2]+")**2")
eq = str(reg.intercept_ - minimum.fun)
for i in range(len(termplot)):
    eq = eq + "+" + str(reg.coef_[i])+"*"+termplot[i]

out = [['Symmetry', symmetry],
       ['Expansion', Doring[symmetry]],
       ['Coefficients', reg.coef_],
       ['R^2', reg.score(terms,energy[:,2])],
       ['Easy Axis', easy],
       ['Hard Axis', hard],
       ['Maximum Energy', -1*maximum.fun - minimum.fun],
       ]

with open('output', 'w') as f:
    f.write(tabulate(out))
    f.write('\n\nPolar Equation:\n{}'.format(eq))
