#!/usr/bin/env python3

###!/usr/bin/env nix-shell
###!nix-shell -i python -p "python37.withPackages(ps: with ps; [ numpy toolz matplotlib])"

import math, re, optparse, operator, os, glob
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call


def load(file):
    f = open(file, 'r')
    r = []
    a = []
    d = []
    for l in f.readlines():
        try:
            n = list(map(str,l.replace(',',' ').split()))
            
            if len(n)>0:
                if n[0].startswith('R'):
                    r.append(n)
                elif n[0].startswith('A'):
                    a.append(n)
                elif n[0].startswith('D'):
                    d.append(n)
        except ValueError:
            pass
    return (np.array(r), np.array(a), np.array(d)) 



def readzmat(filename):
    zmatf = open(filename, 'r')
    atomnames = []
    rconnect = []
    rlist = []
    aconnect = []
    alist = []
    dconnect = []
    dlist = []
    variables = {}
    
    if not zmatf.closed:
        for line in zmatf:
            words = line.split()
            eqwords = line.split('=')
            
            if len(eqwords) > 1:
                varname = str(eqwords[0]).strip()
                try:
                    varval  = float(eqwords[1])
                    variables[varname] = varval
                except:
                    print("Invalid variable definition: " + line)
            
            else:
                if len(words) > 0:
                    atomnames.append(words[0])
                if len(words) > 1:
                    rconnect.append(int(words[1]))
                if len(words) > 2:
                    rlist.append(words[2])
                if len(words) > 3:
                    aconnect.append(int(words[3]))
                if len(words) > 4:
                    alist.append(words[4])
                if len(words) > 5:
                    dconnect.append(int(words[5]))
                if len(words) > 6:
                    dlist.append(words[6])
    
    #replace_vars(rlist, variables)
    #replace_vars(alist, variables)
    #replace_vars(dlist, variables)
    
    return (atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist) 


def write_xyz(atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist):
    npart = len(atomnames)
    
    #print(npart)
    #print('INSERT TITLE CARD HERE')

    xyzarr = np.zeros([npart, 3])
    if (npart > 1):
        xyzarr[1] = [rlist[0], 0.0, 0.0]

    if (npart > 2):
        i = rconnect[1] - 1
        j = aconnect[0] - 1
        r = rlist[1]
        theta = alist[0] * np.pi / 180.0
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        a_i = xyzarr[i]
        b_ij = xyzarr[j] - xyzarr[i]
        if (b_ij[0] < 0):
            x = a_i[0] - x
            y = a_i[1] - y
        else:
            x = a_i[0] + x
            y = a_i[1] + y
        xyzarr[2] = [x, y, 0.0]

    for n in range(3, npart):
        r = rlist[n-1]
        theta = alist[n-2] * np.pi / 180.0
        phi = dlist[n-3] * np.pi / 180.0
        
        sinTheta = np.sin(theta)
        cosTheta = np.cos(theta)
        sinPhi = np.sin(phi)
        cosPhi = np.cos(phi)

        x = r * cosTheta
        y = r * cosPhi * sinTheta
        z = r * sinPhi * sinTheta
        
        i = rconnect[n-1] - 1
        j = aconnect[n-2] - 1
        k = dconnect[n-3] - 1
        a = xyzarr[k]
        b = xyzarr[j]
        c = xyzarr[i]
        
        ab = b - a
        bc = c - b
        bc = bc / np.linalg.norm(bc)
        nv = np.cross(ab, bc)
        nv = nv / np.linalg.norm(nv)
        ncbc = np.cross(nv, bc)
        
        new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z
        new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z
        new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z
        xyzarr[n] = [new_x, new_y, new_z]
            
    #for i in range(npart):
    #    print('{:<4s}\t{:>11.5f}\t{:>11.5f}\t{:>11.5f}'.format(atomnames[i], xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))

    return xyzarr

def set_coordinates(rinlist,ainlist,dinlist, rlist, alist, dlist):
    r = []
    a = []
    d = []


    if len(rinlist) == 0:
        r = [0]*len(rlist)
    else:
        for i in range(len(rlist)):
            key = False
            for j in range(len(rinlist)):
                if rlist[i] == rinlist[j][0]:
                    r.append(float(rinlist[j][1]))
                    key = True
                if j == (len(rinlist)-1) and key == False:
                    r.append(0)

    if len(ainlist) == 0:
        a = [0]*len(alist)
    else:
        for i in range(len(alist)):
            key = False
            for j in range(len(ainlist)):
                if alist[i] == ainlist[j][0]:
                    a.append(float(ainlist[j][1]))
                    key = True
                if j == (len(ainlist)-1) and key == False:
                    a.append(0.0)

    if len(dinlist) == 0:
        a = [0]*len(dlist)
    else:
        for i in range(len(dlist)):
            key = False
            for j in range(len(dinlist)):
                if dlist[i] == dinlist[j][0]:
                    d.append(float(dinlist[j][1]))
                    key = True
                if j == (len(dinlist)-1) and key == False:
                    d.append(0.0)

    return (r, a, d)

def generate_coord(Coord,NM,Q):
    NewCoord = []
    for j in range(len(Coord)):
        NewCoord.append(Coord[j] + (NM[j]*Q))
        
    return NewCoord








    
def step(NM, N, Qi, Qf): # Divide the Normal Coordinates by the Number of Steps
    Q=Qf-Qi
    new=[0]*len(NM)
    for i in range(len(NM)):
        new[i] = NM[i]*Q /( N-1 )
    return new



# Define the Initial Coordinates
def init_coords(type, input, NM, Q):
    Coord=[0]*len(input)

    for i in range(len(input)):
        Coord[i] = input[i] + (NM[i]*Q)
           
    return Coord

# Set the Coordinate for i Step
def step_coords(Coord, NM, i):
    NewCoord = [] 
    for j in range(len(Coord)):
        NewCoord.append(Coord[j] + (NM[j]*i))
    
    return NewCoord


# Set the Input Name
def InputName(file, i):
    if i > 99:
        var = "-%s.input" % i
        input= file.replace(".coord", var);
    elif i > 9:
        var = "-0%s.input" % i
        input= file.replace(".coord", var);            
    else:
        var = "-00%s.input" % i
        input= file.replace(".coord", var);     
    return input



# Write Input File
def WriteInput(input, bs, symb, NewCoord, sym):
    # Write in Input File
    inp = open(input, 'w')
    if sym:
        inp.write("&GATEWAY \nExpert \nSymmetry \n %s \nBasis set \n" % sym)
    else:
        inp.write("&GATEWAY \nExpert \nBasis set\n")
    h=1
    c=1
    o=1
    n=1
    f=1
    b=1
    for j in range(len(NewCoord)):
        # End of Basis
        if j != 0 and symb[j] != symb[j-1]:
            inp.write("End of basis\nBasis set\n")
                
        if symb[j]=="C":
            if c==1:
                inp.write("C.%s\n" % bs)
                inp.write(" %s1    %s Angstrom \n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))
                c=2
            else:
                inp.write(" %s%s    %s Angstrom \n" % (''.join(map(str,symb[j])), c, '  '.join(map(str,NewCoord[j]))))
                c=c+1
        elif symb[j]=="O":
            if o==1:
                inp.write("O.%s\n" % bs)
                inp.write(" %s1    %s Angstrom\n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))
                o=2
            else:
                inp.write(" %s%s    %s Angstrom\n" % (''.join(map(str,symb[j])), o, '  '.join(map(str,NewCoord[j]))))
                o=o+1
        elif symb[j]=="H":
            if h==1:
                inp.write("H.%s\n" % bs)
                inp.write(" %s1    %s Angstrom\n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))
                h=2
            else:
                inp.write(" %s%s    %s Angstrom\n" % (''.join(map(str,symb[j])), h, '  '.join(map(str,NewCoord[j]))))
                h=h+1
        elif symb[j]=="N":
            if n==1:
                inp.write("N.%s\n" % bs)
                inp.write(" %s1    %s Angstrom\n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))
                n=2
            else:
                inp.write(" %s%s    %s Angstrom\n" % (''.join(map(str,symb[j])), n, '  '.join(map(str,NewCoord[j]))))
                n=n+1
        elif symb[j]=="F":
            if f==1:
                inp.write("F.%s\n" % bs)
                inp.write(" %s1    %s Angstrom\n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))
                f=2
            else:
                inp.write(" %s%s    %s Angstrom\n" % (''.join(map(str,symb[j])), f, '  '.join(map(str,NewCoord[j]))))
                f=f+1
        elif symb[j]=="B":
            if b==1:
                inp.write("B.%s\n" % bs)
                inp.write(" %s1    %s Angstrom\n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))
                b=2
            else:
                inp.write(" %s%s    %s Angstrom\n" % (''.join(map(str,symb[j])), b, '  '.join(map(str,NewCoord[j]))))
                b=b+1
        elif symb[j]=="X":
            with open("rydberg.basis") as fmolcas:
                inp.write(fmolcas.read())
                inp.write(" %s    %s Angstrom\n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))
                fmolcas.close()
        else:
            print ("Unknown %s Element!!!" % symb[j])
            exit(0)
                
        # Closing End of Basis
        if j == (len(NewCoord)-1):
            inp.write("End of basis\n")

    inp.write("&SEWARD &END \n NoDelete \nEnd of input\n")
    with open("molcas.input") as fmolcas:
        inp.write(fmolcas.read())
        fmolcas.close()
        inp.close()
            

# Write Molecular Movie
def WriteMovie(inp, symb, NewCoord):
    inp.write(" %s \n\n" % len(NewCoord) )
    for j in range(len(NewCoord)):
        inp.write(" %s    %s \n" % (''.join(map(str,symb[j])), '  '.join(map(str,NewCoord[j]))))



def get_energy(level):
    E=[]
    c=0
    files=sorted(glob.iglob('*.log')) # Used to get all files in numerical order
    for file in files:
        energy=0.0
        for i in open( file ).readlines():
            
            if level == "RASSCF":
                if re.search(r"::    RASSCF root number", i) is not None: # Find energy in .log
                    words = i.split()
                    energy = float( words[7] )  # Energy is the sixth word
                    E.append(energy)
                    
            elif level == "CASPT2":
                if re.search(r"::    CASPT2 Root", i) is not None: # Find energy in .log
                    words = i.split()
                    energy = float( words[6] )  # Energy is the sixth word
                    E.append(energy)
                    
            elif level == "MS-CASPT2":
                if re.search(r":    MS-CASPT2 Root", i) is not None: # Find energy in .log
                    words = i.split()
                    energy = float( words[6] )  # Energy is the sixth word
                    E.append(energy)
                
            elif level == "RASSI":
                if re.search(r"::    RASSI State ", i) is not None: # Find energy in .log
                    words = i.split()
                    energy = float( words[6] )  # Energy is the sixth word
                    E.append(energy)
            else:
                print("You forgot something, right? Maybe... level of calculation???")
                exit()
            
        if energy == 0.0:
            E.append(energy)
            
        c=c+1
        if c == 1:
            nstates=len(E)
    return E, nstates


# MAIN PROGRAM
def main():
    import sys
    f = optparse.OptionParser()
    # Get Cartesian Coordinates
    f.add_option('-z', '--zmat', type = str, default = None)
    # Get Normal Coordinates
    f.add_option('-n', '--nm' , type = str, default = 'co2-nm.coord')
    # Get Normal Coordinates
    f.add_option('-e', '--eq' , type = str, default = None)    
    # Get Number of Steps
    f.add_option('-N', '--N'  , type = int, default = 10)
    # Get Initial Amplitude
    f.add_option('-i', '--qi'  , type = float, default = -3)
    # Get Final Amplitude
    f.add_option('-f', '--qf'  , type = float, default = 3)
    # Get Type of Run
    f.add_option('-t', '--typ', type = str, default = 'help')
    # Get Number of States
    f.add_option('-s', '--st' , type = int, default = 1)    
    # Get Number of Atom from Amplitude Calculation (Only for typ=get)
    f.add_option('-a', '--atm' , type = int, default = 2)
    # Basis set
    f.add_option( '-b', '--bs' , type = str, default = 'ANO-RCC-VTZP')
    # Symmetry
    f.add_option( '-c', '--sym' , type = str, default = None)
    (arg, args) = f.parse_args(sys.argv[1:])




    if arg.typ == "pes":

        atomnames,rconnect,rlist,aconnect,alist,dconnect,dlist=readzmat(arg.zmat)
        reqlist,aeqlist,deqlist=load(arg.eq)
        rnmlist, anmlist, dnmlist=load(arg.nm)
        
        req,aeq,deq=set_coordinates(reqlist,aeqlist,deqlist, rlist, alist, dlist)
        rnm,anm,dnm=set_coordinates(rnmlist,anmlist,dnmlist, rlist, alist, dlist)

        
        # Get Step Ratio
        if (arg.N % 2 == 0):    # The N. Steps must be even, to have the Equilibrium
            arg.N = arg.N + 1   # Geometry in the center of the Curve
        rstp = step(rnm, arg.N, arg.qi, arg.qf)
        astp = step(anm, arg.N, arg.qi, arg.qf)
        dstp = step(dnm, arg.N, arg.qi, arg.qf)

        # Get Initial Coordinates
        
        InitR = init_coords("R",req, rstp, arg.qi)
        InitA = init_coords("A",aeq, astp, arg.qi)
        InitD = init_coords("D",deq, dstp, arg.qi)

        
        # Write XYZ Files
        file=arg.nm
        inp = open("movie.dat", 'w')
        rsize=arg.N

        Q=np.linspace(arg.qi,arg.qf,arg.N)
        
        for i in range(arg.N):
            
            NewR=generate_coord(req, rnm, Q[i])
            NewA=generate_coord(aeq, anm, Q[i])
            NewD=generate_coord(deq, dnm, Q[i])           
            
           
            input=InputName(file,i)
        
            xyz=write_xyz(atomnames, rconnect, NewR, aconnect, NewA, dconnect, NewD)

    
            WriteInput(input, arg.bs, atomnames, xyz, arg.sym)
            # Write Movie
            WriteMovie(inp, atomnames, xyz)

            job=input.replace(".input", ".sh")

            scpt=open(job, 'w')

            with open("submit.sh") as fugu:
                scpt.write(fugu.read())
            scpt.write("pymolcas -b 1 -f %s" % input)
            fugu.close()
            scpt.close()

            #call("sbatch %s" % job, shell=True)




            
    # GET RESULTS
    elif arg.typ == "get":
        print ("----------------------------------- GET PEC MODULE -------------------------------\n")
        
        # Get energy from .log Files
        for i in range(arg.st):
            energy=PEC.get_energy(i+1)

        # Get Amplitude Q from .xyz Files
        amp = PEC.get_amplitude(arg.xyz, arg.nm, arg.atm)
        
        print ("      Results\n Amp.   Energy(au)")
        x=[]
        y=[]
        for i in range(len(amp)):
            if energy[i]==0.0:
                pass
            else:
                x.append(amp[i])
                y.append(energy[i])

        for i in range(len(x)):
            print (x[i], y[i])

        if  arg.st == 1:
            plt.title('Ground State', fontsize=22)
        else:
            plt.title('Excited State', fontsize=22)
        plt.xlabel('Amplitude', fontsize=18)
        plt.ylabel('Energy (a.u.)', fontsize=18)
        plt.plot(x,y)
        plt.show()

    else:
        print ("\n WRONG type! Please, check help information (-t help)!\n")

    print ("\n DONE!\n Have a Nice Day!\n")

if __name__=="__main__":
    main()
