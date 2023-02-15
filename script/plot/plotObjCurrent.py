
import h5py
import pylab as plt
import numpy as np

k_B=1.380649*10**(-23) # boltzman constant SI units

Geometry= "Sphere"#"None" #"Cylinder"

#If using a geometry define params for OML current
#electron current
n=5.8977e9 #m^-3
q=-1.60217662e-19 #C
m= 9.10938356e-31 #Kg
V_th=123111 #m/s
T=(m*V_th*V_th)/k_B #m/s
V=4.5 #V

#probe shape and current
#stepsize=0.0035
stepsize=0.007 #0.0035*2
radius_N_cells=2 # in meters radius = radius_N_cells*stepsize

#For Cylinder
length=stepsize*59

def calc_EMA(array,S_val,N_val):

    # S_val: Smothing factor
    # N_val: number of point to use in filter

    new_array=np.zeros(len(array))
    k=S_val/(N_val+1)
    new_array[0]=array[0]
    for t in range(1,len(array)):
        new_array[t]=array[t]*k+new_array[t-1]*(1-k)

    return new_array


def calc_OML(n,q,S,T,m,V,K,beta):
	I_th=n*q*S*np.sqrt((k_B*T)/(2*np.pi*m))
	eta=-(q*V)/(k_B*T)
	print('dimensionless probe potentioal eta = %.1f'%eta)
	I=I_th*K*(1+eta)**beta
	return I

def calc_OML_sphere(r,n,q,T,m,V,arr_len):
	S=4*np.pi*r*r
	K=1
	beta=1
	I=np.zeros(arr_len)
	I.fill(calc_OML(n,q,S,T,m,V,K,beta))
	print(q)
	print(calc_OML(n,q,S,T,m,V,K,beta))
	return I

def calc_OML_cylinder(h,r,n,q,T,m,V,arr_len):
        S=2*np.pi*r*(h+r)
        K=2/np.sqrt(np.pi)
        beta=0.5 # 0.5 for zero end effects
        I=np.zeros(arr_len)
        I.fill(calc_OML(n,q,S,T,m,V,K,beta))
        print(q)
        print(calc_OML(n,q,S,T,m,V,K,beta))
        return I

dataset='data'
hist = h5py.File('../../'+dataset+'/history.xy.h5','r')
curr1 = hist['/current/electrons/dataset']
curr2 = hist['/current/ions/dataset']

#hist2 = h5py.File('../../data_nocoll/history.xy.h5','r')
#curr3 = hist2['/current/electrons/dataset']
#curr4 = hist2['/current/ions/dataset']
#print(curr1)

curr1 = curr1[:,1]		# Extract y-axis
curr2 = curr2[:,1];
#curr3 = curr3[:,1]		# Extract y-axis
#curr4 = curr4[:,1];

# Uncoment these to use Exponential Moving Average
curr1=calc_EMA(curr1,2,100)
curr2=calc_EMA(curr2,2,100)



if Geometry=="Sphere":
	#Spehere
	radius=stepsize*radius_N_cells #most probable
	curr_OML=calc_OML_sphere(radius,n,q,T,m,V,len(curr1))
	radius=stepsize*np.sqrt((radius_N_cells-1)**2 +(radius_N_cells-1)**2 ) #min
	curr_OML_min=calc_OML_sphere(radius,n,q,T,m,V,len(curr1))
	radius=stepsize*np.sqrt((radius_N_cells)**2 +(radius_N_cells-1)**2 ) #max
	curr_OML_max=calc_OML_sphere(radius,n,q,T,m,V,len(curr1))
	plt.plot(curr_OML,label='OML current most likely')
	plt.plot(curr_OML,label='OML current max')
	plt.plot(curr_OML_min,label='OML current min')

	print("error last step is %.1f "%(((curr1[-1]-curr_OML[-1])/curr_OML[-1])*100) )


if Geometry=="Cylinder":
	#Cylinder

	radius=stepsize*radius_N_cells #most probable
	curr_OML=calc_OML_cylinder(length,radius,n,q,T,m,V,len(curr1))
	radius=stepsize*np.sqrt((radius_N_cells-1)**2 +(radius_N_cells-1)**2 ) #min
	curr_OML_min=calc_OML_cylinder(length,radius,n,q,T,m,V,len(curr1))
	radius=stepsize*np.sqrt((radius_N_cells)**2 +(radius_N_cells-1)**2 ) #max
	curr_OML_max=calc_OML_cylinder(length,radius,n,q,T,m,V,len(curr1)) #max
	
	plt.plot(curr_OML,label='OML current most likely')
	plt.plot(curr_OML,label='OML current max')
	plt.plot(curr_OML_min,label='OML current min')

	print("error last step is %.1f "%(((curr1[-1]-curr_OML[-1])/curr_OML[-1])*100) )


plt.plot(curr1,label='electron current')
plt.plot(curr2,label='ion current')


plt.legend(loc='lower right')
plt.ylabel("I [A]")
plt.xlabel("timesteps")
#plt.savefig("data/currents_both.png")
plt.show()

