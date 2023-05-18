
from fnmatch import fnmatch, fnmatchcase

def get_order(filename):
    order=int(open(filename.split(".")[0]).readlines()[0].split("\t")[0])
    return order



def get_Li_added(filename):
    
    atomorder=get_order(filename)
    ##need specify filename in mol and atomorder.

    atom_order=atomorder-1 ####(原子序号减去1)
    
    f=open(filename).readlines()
    row13=f[0:3]
    an=int(f[3].split()[0])
    
    
    row4=str(an+1)+"  "+f[3][4:-1]+"\n"
    
    rowleft=f[4+an:-1]  
    
    #Lix=
    #Liy=
    #Liz=
    
    #rowLi='    6.7123   -3.0161    0.2849 Li   0  0  0  0  0  0  0  0  0  0  0  0\n'
    
    
    atomrow=f[4:4+an]
    
    
    
    
    
    center=list(map(float,atomrow[atom_order].split()[0:3]))
    
    dist_list=[10]
    
    for i in range(an):
        atom_coord=list(map(float,atomrow[i].split()[0:3]))
        
        x2=atom_coord[0]-center[0]
        y2=atom_coord[1]-center[1]
        z2=atom_coord[2]-center[2]
        
        dist=(x2**2+y2**2+z2**2)**0.5
        
        if dist < min(dist_list) and dist!=0:
            dist_list.append(dist)
            min_index=i
       
    nearest_atom=list(map(float,atomrow[min_index].split()[0:3]))
    
    a=(1/min(dist_list))*2
    lix=-(nearest_atom[0]-center[0])*a+center[0]
    liy=-(nearest_atom[1]-center[1])*a+center[1]
    liz=-(nearest_atom[2]-center[2])*a+center[2]
    
    #print(lix,liy,liz)
    
    
    
    rowLi='    '+str(round(lix,4))+'   '+str(round(liy,4))+'    '+str(round(liz,4))+' Li   0  0  0  0  0  0  0  0  0  0  0  0\n'
    
    
    
    atomrow.append(rowLi)
    
    #print(rowLi)
    
    with open("Li_added"+filename,'w+') as f:
        f.writelines(row13)
        f.writelines(row4)
        f.writelines(atomrow)
        f.writelines(rowleft)





get_Li_added("mol_88.mol")


