import numpy as np
import functions,os

def gen_cfg(f_min,f_max,f_step,alpha_min,alpha_max,alpha_step):

    opath = "cfg_folder"

    if not os.path.exists(opath):
        os.system("mkdir "+opath)
    os.system("rm %s/*"%opath)

    faxis = np.arange(f_min,f_max+f_step/2,f_step)
    alpha_axis = np.arange(alpha_min,alpha_max+alpha_step/2,alpha_step)

    i = 0
    for f in faxis:
        for alpha in alpha_axis:
            o = open("%s/cfg_%d"%(opath,i),"w")
            functions.write_table([[f],[alpha]],o)
            o.close()
            i += 1

if __name__ == "__main__":
    gen_cfg(0,0.1,0.05,0,90,45)
            
