

filename="mol_998-pos-1.xyz"


def get_last_frame(filename):
    f=open(filename).readlines()
    
    for i in range(len(f)):
        
        if "=" in f[i]:
            index=i
    
    
    text=f[index-1:]
    
    print(text)
    with open("Opted-"+filename,"w+") as f:
        f.writelines(text)
            
        
get_last_frame(filename)