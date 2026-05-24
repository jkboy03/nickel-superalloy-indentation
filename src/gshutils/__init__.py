# from functools import cache 
import numpy as np
from scipy.optimize import linprog


class gshUtils():
    def __init__(self,symmetry,truncation):
        if symmetry == "cubic_triclinic":
            from .cubic_triclinic_wrapped_tfunction import governor
        else:
            print(f"{symmetry} is not a supported symmetry")
            
        self.funcs, self.indx = governor()
        
        self.nStates = np.sum(self.indx[:,0] <= truncation)
            
    #convert one angle to gsh
#     @cache
    def angle_to_gsh(self,angle0,angle1,angle2):
        return np.array([self.funcs[i](angle0,angle1,angle2) for i in range(0,self.nStates)])

    #split the angles up from an array to individual values. allows the cache to work. make fast
    def split_angs_hash(self,angles):
        return self.angle_to_gsh(angles[0],angles[1],angles[2])
    
    #convert one angle to ssh
#     @cache
    def angle_to_ssh(self,angle0,angle1,angle2):
        return np.array([(1/np.sqrt((4*np.pi)/(2*self.indx[i][0]+1)))*self.funcs[i](angle0,angle1,angle2)
                         for i in range(0,self.nStates) if self.indx[i][1]==0])
    
    def split_angs_hash_ssh(self,angles):
        return self.angle_to_ssh(0,angles[0],angles[1]+np.pi/2)

    #convert a struct of euler angles to a struct of gsh coeffs
    def struct_to_gsh(self,eulerAngles):
        return np.apply_along_axis(self.split_angs_hash, -1, eulerAngles)
    
    def struct_to_ssh(self,eulerAngles):
        return np.apply_along_axis(self.split_angs_hash_ssh, -1, eulerAngles)

    #sorry to anyone trying to decipher the next two functions. It's either write 500 line of code or make them these nuanced rascals. you just gotta live with it...
    def gsh_complex_to_gsh_mixed(self,gshComplex):
        gshMixed = np.zeros_like(gshComplex)

        for i in range(self.nStates):
            if self.indx[i,1]==0:
                gshMixed[...,i]=gshComplex[...,i].real
            elif self.indx[i,1]<0:
                gshMixed[...,i]=gshComplex[...,i].imag
            elif self.indx[i,1]>0:
                gshMixed[...,i]=gshComplex[...,i].real

        return gshMixed.real

    #warning wont work if you dont truncate at changes in L
    def gsh_mixed_to_gsh_complex(self,gshMixed):
        gshComplex = np.zeros_like(gshMixed).astype("complex")

        for i in range(self.nStates):
            if self.indx[i,1]==0:
                gshComplex[...,i] = gshMixed[...,i] + 0j
            elif (self.indx[i,1] % 2 != 0) and self.indx[i,1]<0:
                gshComplex[...,i] = -1*gshMixed[...,i-2*self.indx[i,1]] + 1j*gshMixed[...,i]

            elif (self.indx[i,1] % 2 == 0) and self.indx[i,1]>0:
                gshComplex[...,i] = gshMixed[...,i] - 1j*gshMixed[...,i-2*self.indx[i,1]]
            elif self.indx[i,1]<0:
                gshComplex[...,i]+= 1j*gshMixed[...,i] + gshMixed[...,i-2*self.indx[i,1]]
            elif self.indx[i,1]>0:
                gshComplex[...,i] = gshMixed[...,i] + 1j*gshMixed[...,i-2*self.indx[i,1]]

        return gshComplex
    
    def __p2_crosscorrelation(self,arr1, arr2):
        ax = list(range(0, len(arr1.shape)))
        arr1_FFT = np.fft.rfftn(arr1, axes=ax)
        arr2_FFT = np.fft.rfftn(arr2, axes=ax)
        return np.fft.irfftn(arr1_FFT.conjugate() * arr2_FFT, s=arr1.shape, axes=ax).real / np.product(
            arr1.shape)

    def gsh_statistics(self,struct):
        stats = np.zeros_like(struct, dtype=np.float64)
        for i in range(0, struct.shape[-1]):
            stats[..., i] = self.__p2_crosscorrelation(struct[..., 1], struct[..., i])
        return (stats[...,1:], stats[...,0])
    
    def volume_weighted_mean(self,struct):
        return np.mean(struct.reshape(-1,self.nStates),axis=0)
    
    def plot_odf(self,odf_coeffs):
        plot_angles = np.ones((50,3)) ##todo 
        basis = np.array([self.funcs[i](plot_angles[:,0],plot_angles[:,1],plot_angles[:,2]) for i in range(0,self.nStates)])
        value = np.dot(odf_coeffs,basis)
        
        
        return 
    
    
    
    
    
    #################idk if these are useful anymore
    
    def __single_compare(self,P1,P2,target):

        if np.sum(P2-P1) == 0:
            return (2e300,None,None)

        t = np.sum((P2-P1)*(target-P1))/np.sum((P2-P1)**2)
        nearest_point = P1+(P2-P1)*t
        distance = np.linalg.norm(nearest_point-target)
        weights = (1-t,t)

        return (distance,nearest_point,weights)

    def __best_pair(self,Start,Positions,target,children_list,enforce_children=False):
        distance_min = 1e300
        for i,Position in enumerate(Positions):
            results = self.__single_compare(Start,Position,target)
            if (results[0] < distance_min) and (results[2][1]>0) and (results[2][1]<1):
                    distance_min = results[0]
                    best_results = results
                    best_pos = i


        return (best_pos,*best_results)

    def weight_convex_set(self,Positions,target,threshold,max_iter):

        weights_list = np.zeros(len(Positions))

        children_list=[]

        start_index  = 0 

        pair_index,distance,nearest_point,weights = self.__best_pair(Positions[start_index],Positions,target,children_list)

        weights_list[start_index] =  weights[0]
        weights_list[pair_index] =  weights[1]
        children_list = [start_index,pair_index]

        for i in range(max_iter):
            pair_index,distance,nearest_point,weights = self.__best_pair(nearest_point,Positions,target,children_list)
            for child in children_list:
                weights_list[child] *= weights[0]

            if weights_list[pair_index] == 0:
                weights_list[pair_index] =  weights[1]
            else:
                weights_list[pair_index] += weights[1]

            children_list.append(pair_index)

            children_list = list(set(children_list))

            if distance < threshold:
                print(f'[Info] Converged in {i} Iterations')
                return weights_list
            
        print("[Warning] Did not Converge...")
        return weights_list
    
    def weight_convex_set_scipy(self,Positions,target,tol):
        linprog(Positions,
                method='simplex',
                options=
                {'maxiter': 5000, 
                 'disp': True, 
                 'presolve': True, 
                 'tol': tol, 
                 'autoscale': False, 
                 'rr': True, 
                 'bland': False},
                x0=None)
    


