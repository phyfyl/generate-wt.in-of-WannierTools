#这个脚本会自动生成需要计算的任务的wanniertools输入文件的模板以及提交任务的脚本
#运行这个脚本会在主文件夹下生成一个"wt_calc"的文件夹，所有wanniertools的计算任务的文件夹都会在这个文件下，
#不同的任务会有不同的文件夹
#输入参数
queue = "e52680v3ib!"    #队列名
num_cores = 48           #使用的核数
job_name = ""            #任务名,如果有输入就用输入的，无输入默认用“主文件夹名_任务名”的形式
wt_path = "wt_fyl.x"     #wt.x路径
soc=True                 #是否考虑SOC, True or Flase
calc_task=24             #计算任务，1-24
#calc_task值对应的任务：
# 1:Bulk band calculation
# 2:3D Fermi Surface calculation
# 3:A cross-section of the Fermi surface
# 4:Bulk spin texture calculation
# 5:Density of states calculation
# 6:Find Nodes calculation
# 7:Energy gap calculations_plane
# 8:Energy gap calculations_cube
# 9:Slab band calculation
# 10:Nanowire/nanoribbon band calculation
# 11:Surface states ARPES calculation
# 12:Surface states QPI calculation
# 13:Fermi arc calculation
# 14:Berry phase calculation
# 15:Berry curvature calculation for 3D bulk case
# 16:Anomalous Hall conductivity calculation
# 17:Wilson loop calculation
# 18:Z2 number for 3D bulk case
# 19:Chern number for 3D bulk case
# 20:Mirror Chern number calculation
# 21:Weyl Chirality calculation
# 22:Landau level calculations
# 23:Ordinary magnetoresistance calculations
# 24:CPGE caculation

import os
import numpy as np

task_labels = {1:'BulkBand_calc = T',2:'BulkFS_calc = T',3:'BulkFS_plane_calc = T\nTranslate_to_WS_calc  = F',4:'BulkSpintexture_calc = T',\
               5:'DOS_calc = T',6:'FindNodes_calc = T',7:'BulkGap_Plane_calc = T',8:'BulkGap_Cube_calc = T',\
               9:'SlabBand_calc = T',10:'WireBand_calc = T',11:'SlabSS_calc = T',12:'SlabQPI_kplane_calc = T',\
               13:'SlabArc_calc=T',14:'BerryPhase_calc = T',15:'BerryCurvature_calc=T',16:'AHC_calc=T',\
               17:'WannierCenter_calc=T',18:'Z2_3D_calc = T',19:'Chern_3D_calc = T',20:'MirrorChern_calc=T',\
               21:'WeylChirality_calc = T',22:'BulkBand_calc = T\nHof_Butt_calc = T\nLandauLevel_k_calc = T\nLandauLevel_wavefunction_calc = F',\
               23:'Boltz_OHE_calc= T\nSymmetry_Import_calc  = T',24:'cpge_calc = T'}
task_names = {1:'BulkBand_calc',2:'BulkFS_calc',3:'BulkFS_plane_calc',4:'BulkSpintexture_calc',\
               5:'DOS_calc',6:'FindNodes_calc',7:'BulkGap_Plane_calc',8:'BulkGap_Cube_calc',\
               9:'SlabBand_calc',10:'WireBand_calc',11:'SlabSS_calc',12:'SlabQPI_kplane_calc',\
               13:'SlabArc_calc',14:'BerryPhase_calc',15:'BerryCurvature_calc',16:'AHC_calc',\
               17:'WannierCenter_calc',18:'Z2_3D_calc',19:'Chern_3D_calc',20:'MirrorChern_calc',\
               21:'WeylChirality_calc',22:'LandauLevel_calc',23:'Boltz_OHE_calc',24:'cpge_calc'}

kcube_bulk_label = [2,5,6,8,16,23,24]
kplane_bulk_label = [3,4,7,15,17,20]
kpath_bulk_label = [1,22]
kpath_slab_label = [9,11]
kplane_slab_label = [12,13]
Nk1_label = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
Nk2_label = [2,3,4,5,6,7,12,13,15,16,17,18,19,20,21,23,24]
Nk3_label = [2,5,6,8,16,23,24]
Eta_arc_label = [3,4,5,11,12,13,24]
E_arc_label = [3,4,22]
SELECTED_WANNIERORBITALS_label = [1]
SELECTED_ATOMS_label = [4]
Omega_label = [5,11,16,22,23,24]
Gap_threshold_label = [6]
NSLAB_label = [9,22]
NSLAB1_label = [10]
NSLAB2_label = [10]
NP_label = [11]
KPATH_BERRY_label = [14]
WEYL_CHIRALITY_label = [21]
SURFACE_label = [9,10,11,12,13,22]
B_direction_label = [23]
B_tau_label = [23]
T_label = [23]
SELECTEDORBITALS_label = [22]
SELECTEDBANDS_label = [23]
freq_label = [24]
delta_op_label = [24]

def read_poscar(poscar_file):
    """读取POSCAR文件并提取晶格信息和原子坐标"""
    with open(poscar_file, 'r') as f:
        lines = [line.strip() for line in f.readlines()]

    # 解析晶格常数和基矢
    scale = float(lines[1])
    lattice = []
    for i in range(2, 5):
        lattice.append([float(x)*scale for x in lines[i].split()[:3]])
    
    # 解析原子信息
    elements = lines[5].split()
    num_atoms = list(map(int, lines[6].split()))
    sum_atoms = sum(num_atoms)
    coord_type = lines[7].lower()
    
    # 解析原子坐标
    atoms = []
    index = 8
    for elem, num in zip(elements, num_atoms):
        for _ in range(num):
            coords = list(map(float, lines[index].split()[:3]))
            index += 1
          #  if coord_type == 'direct':
                # 转换为笛卡尔坐标
          #      x, y, z = [
          #          coords[0]*lattice[0][0] + coords[1]*lattice[1][0] + coords[2]*lattice[2][0],
           #         coords[0]*lattice[0][1] + coords[1]*lattice[1][1] + coords[2]*lattice[2][1],
           #         coords[0]*lattice[0][2] + coords[1]*lattice[1][2] + coords[2]*lattice[2][2]
            #    ]
           # else:  # Cartesian坐标直接缩放
            x, y, z = [scale * c for c in coords]
            atoms.append({"element": elem, "x": x, "y": y, "z": z})
    
    return lattice, num_atoms,atoms

def read_wannier_win(win_file):
    """读取wannier90.win文件并提取kpath和投影轨道信息"""
    with open(win_file, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    # 提取k点路径
  #  kpath, in_kpath = [], False
  #  for line in lines:
  #      if 'begin kpoint_path' in line.lower():
  #          in_kpath = True
  #      elif 'end kpoint_path' in line.lower():
  #          in_kpath = False
  #      elif in_kpath and line:
  #          kpath.append(line)
    
    # 提取投影轨道
    projections, in_proj = [], False
    for line in lines:
        if 'begin projections' in line.lower():
            in_proj = True
        elif 'end projections' in line.lower():
            in_proj = False
        elif in_proj and line and ':' in line:
            coord_part, orb_part = line.split(':', 1)
          #  coords = coord_part.strip().split()
            projs=[]
            for orb in orb_part.split(';'):
                orb_clean = orb.strip()
                if orb_clean:
                    projs.append(orb_clean)
            projections.append(projs)
    
    return  projections

def get_fermi_occupation(eigenval,win_file,fermi_file):
    with open(win_file, 'r') as f1:
        lines_win = [line.strip() for line in f1.readlines()]
    with open(eigenval, 'r') as f2:
        lines_eval = [line.strip() for line in f2.readlines()]
    with open(fermi_file, 'r') as f3:
        lines_fermi = [line.strip() for line in f3.readlines()]
    
    fermi=float(lines_fermi[1])
    for line in lines_win:
        if "dis_win_min" in line.lower():
            win_min = float(line.split("=")[1].strip())
    
    num_ocps=0
    for i,line in enumerate(lines_eval):
        data=line.split()
        en = float(data[2].strip())
        if en > win_min and en < fermi:
        #   ocps=float(data[2].strip()
        #    if ocps>0.9:
            num_ocps+=1
        elif en > fermi:
            break


    return fermi,num_ocps
            
def read_kpath(kpath_file):
    kpath=[]
    with open(kpath_file, 'r') as f:
        lines = f.readlines()
    for i,line in enumerate(lines):
        data = line.split()
        kp=[]
        if i>3 and len(data)>0:
            kp.append(data[-1].strip())
            kp.append(data[0].strip())
            kp.append(data[1].strip())
            kp.append(data[2].strip())
            kpath.append(kp)
    print(kpath)
    return kpath
    
def write_wt_in(lattice,num_atoms, atoms, projections, fermi,num_ocps,kpath,soc,output_file='wt.in'):
    """生成WannierTools输入文件"""
    with open(output_file, 'w') as f:
        f.write("&TB_FILE\n")
        f.write("Hrfile = 'wannier90_hr.dat'\n")
        f.write("Particle='electron'\n/\n")

        f.write("\nLATTICE\n")
        f.write("Angstrom\n")
        for vec in lattice:
            f.write(f"{vec[0]:.9f} {vec[1]:.9f} {vec[2]:.9f}\n")
        
        N_atoms=sum(num_atoms)
        f.write("\nATOM_POSITIONS\n")
        f.write(f"{N_atoms}\n")
        f.write("Direct\n")
        for atom in atoms:
            f.write(f"{atom['element']} {atom['x']:.9f} {atom['y']:.9f} {atom['z']:.9f}\n")
        
        # 投影轨道
        f.write("\nPROJECTORS\n")
        for i in range(len(num_atoms)):
            n = num_atoms[i]
            num_projs = len(projections[i])
            for j in range(n):
                f.write(f"{num_projs} ")
        f.write("\n")
        for i in range(len(num_atoms)):
            n = num_atoms[i]
            num_projs = len(projections[i])
            for j in range(n):
                atom=atoms[sum(num_atoms[0:i])+j]
                f.write(f"{atom['element']} ")
                for k in range(num_projs):
                    f.write(f"{projections[i][k]} ")
                f.write("\n")

        # 控制参数
        f.write("\n&CONTROL\n")
        f.write(f"{task_labels[calc_task]}\n")
        f.write("/\n")

        f.write("\n&SYSTEM\n")
        f.write(f"SOC = {soc}\n")
        f.write(f"Numoccupied = {num_ocps}\n")
        f.write(f"E_FERMI = {fermi}\n")
        if calc_task in NSLAB_label:
            f.write("NSLAB = 1\n")
        if calc_task in NSLAB1_label:
            f.write("NSLAB1 = 1\n")
        if calc_task in NSLAB2_label:
            f.write("NSLAB2 = 1\n")
        if calc_task in B_direction_label:
            f.write("Btheta = 0, Bphi = 90   ! magnetic field direction\n")
        f.write("/\n")

        f.write("\n&PARAMETERS\n")
        if calc_task in Eta_arc_label:
            f.write("Eta_Arc = 0.0259\n")
        if calc_task in E_arc_label:
            f.write("E_Arc = 0.0\n")
        if calc_task in Omega_label:
            f.write(f"OmegaNum = 101\n")
            f.write("OmegaMin = 0\n")
            f.write("OmegaMax = 1\n")
        if calc_task in freq_label:
            f.write("freqnum = 101\n")
            f.write("freqmin = 1\n")
            f.write("freqmax = 2\n")
        if calc_task in delta_op_label:
            f.write("delta_op = 0.02\n")
        if calc_task in Nk1_label:
            f.write("Nk1 = 101\n")
        if calc_task in Nk2_label:
            f.write("Nk2 = 101\n")
        if calc_task in Nk3_label:
            f.write("Nk3 = 101\n")
        if calc_task in NP_label:
            f.write("NP = 2\n")
        if calc_task in Gap_threshold_label:
            f.write("Gap_threshold = 0.001  ! a value to determine which point should be identified as a node\n")
        if calc_task in B_tau_label:
            f.write("BtauNum = 100    ! Number of B*tau we calculate\n")
            f.write("BtauMax = 40.0   ! The maximum B*tau, starting from Btau=0\n")
        if calc_task in T_label:
            f.write("Tmin = 30        ! Temperature in Kelvin\n")
            f.write("Tmax = 330       ! Temperature in Kelvin\n")
            f.write("NumT = 11        ! number temperature we calculate. T_i=Tmin+(Tmax-Tmin)/(NumT-1)*(i-1)\n")
        f.write("/\n")

        if calc_task in SURFACE_label:
            f.write("\nSURFACE\n")
            f.write("1 0 0\n")
            f.write("0 1 0\n")
            f.write("0 0 1\n")

        if calc_task in SELECTED_WANNIERORBITALS_label:
            f.write("\nSELECTED_WANNIERORBITALS  ! get projected bands onto different orbitals, here we only consider orbital and omit the spin freedom\n")
            f.write("2\n")
            f.write("1-6\n")
            f.write("7-15\n")

        if calc_task in SELECTEDORBITALS_label:
            f.write("\nSELECTEDORBITALS\n")
            f.write("1   ! NumberofSelectedOrbitals without spin degeneracy\n")
            f.write("1   ! SelectedOrbitals indices without spin degeneracy\n")

        if calc_task in SELECTEDBANDS_label:
            f.write("\nSELECTEDBANDS\n")
            f.write("1\n")
            f.write("6  ! the 6'th band is crossing the Fermi level.\n")

        if calc_task in SELECTED_ATOMS_label:
            f.write("\nSELECTED_ATOMS  !projection only onto the selected atoms\n")
            f.write("2 ! number groups of selected atoms\n")
            f.write("6 12 18 24 30  ! top surface's atoms\n")
            f.write("1  7 13 19 25  ! bottom surface's atoms\n")

        if calc_task in kcube_bulk_label:
            f.write("\nKCUBE_BULK\n")
            f.write("-0.5 -0.5 -0.5\n")
            f.write("1 0 0\n")
            f.write("0 1 0\n")
            f.write("0 0 1\n")

        if calc_task in kplane_bulk_label:
            f.write("\nKPLANE_BULK\n")
            f.write("0.00 0.00 0.00\n")
            f.write("1.00 0.00 0.00\n")
            f.write("0.00 1.00 0.00\n")

        if calc_task in kpath_bulk_label:
            f.write("\nKPATH_BULK\n")
            n_kpath=round(len(kpath)/2)
            f.write(f"{n_kpath}\n")
            for i in np.arange(0,2*n_kpath,2):
                f.write(f'''{kpath[i][0]:<7} {' '.join([f"{float(ele):>10.6f}" for ele in kpath[i][1:4]])} \
    {kpath[i+1][0]:<7} {' '.join([f"{float(ele):>10.6f}" for ele in kpath[i+1][1:4]])}\n''')

        if calc_task in kpath_slab_label:
            f.write("\nKPATH_SLAB\n")
            f.write("2\n")
            f.write("K 0.33 0.67 G 0.0 0.0\n")
            f.write("G 0.0 0.0 M 0.5 0.5\n")
        
        if calc_task in kplane_slab_label:
            f.write("\nKPLANE_SLAB\n")
            f.write("0.00 0.00 0.00\n")
            f.write("1.00 0.00 0.00\n")
            f.write("0.00 1.00 0.00\n")

        if calc_task in KPATH_BERRY_label:
            f.write("\nKPATH_BERRY\n")
            f.write("11\n")
            f.write("Direct")
            f.write('''
 0.3    0.333  -0.2
 0.3    0.333  -0.1
 0.3    0.333  -0.0
 0.3    0.333   0.1
 0.3    0.333   0.2
 0.33   0.333   0.2
 0.33   0.333   0.1
 0.33   0.333   0.0
 0.33   0.333  -0.1
 0.33   0.333  -0.2
 0.3    0.333  -0.2''')
            
        if calc_task in WEYL_CHIRALITY_label:
            f.write("\nWEYL_CHIRALITY\n")
            f.write("8           ! Num_Weyls\n")
            f.write("Cartesian   ! Direct or Cartesian coordinate\n")
            f.write("0.004       ! Radius of the ball surround a Weyl point")
            f.write('''
 0.219436   -0.045611   -0.000000    ! Positions of Weyl points, No. of lines should larger than Num_weyls
-0.219515   -0.045063   -0.000000
 0.220195   -0.038682   -0.000000
-0.220183   -0.038936   -0.000000
 0.219514    0.045063    0.000000
-0.219434    0.045620    0.000000
-0.220194    0.038678    0.000000
 0.220181    0.038941    0.000000''')

def write_job(path):
    global job_name
    if len(job_name.split())==0:
        job_name=f"{path.split('/')[-1]}_{task_names[calc_task]}"
    with open('job','w') as f:
        f.write(f'''#BSUB -q {queue}
#BSUB -n {num_cores}
#BSUB -J {job_name}
#BSUB -e err
#BSUB -o out

module load ips/2018u4
mpirun {wt_path}''')
        
def main():
    #创建文件夹
    while True:
        dir = os.listdir()
        if "wannier90.win" in dir:
            main_path = os.getcwd()
            if "wt_calc" not in dir:
                os.system("mkdir wt_calc")
                print(f"Created wt_calc directory under {main_path}")
            os.chdir("wt_calc")
            dir_wt=os.listdir()
            if task_names[calc_task] not in dir_wt:
                os.system(f"mkdir {task_names[calc_task]}")
                print(f"Created {task_names[calc_task]} directory under wt_calc")
            os.chdir(f"{task_names[calc_task]}")
            os.system("cp ../../scf/FERMI_ENERGY.in .")
            os.system("cp ../../wannier/POSCAR .")
            os.system("cp ../../wannier/wannier90.win .")
            os.system("cp ../../wannier/wannier90.eig .")
            os.system("cp ../../wannier/wannier90_hr.dat .")
            if calc_task==1 or calc_task==22:
                os.system("cp ../../band/KPOINTS .")
                kpath = read_kpath('KPOINTS')
            else:
                kpath = None
            break
        else:
            os.chdir('..')

    # 检查文件存在性
    if not os.path.exists('POSCAR'):
        raise FileNotFoundError("POSCAR not found in current directory")
    if not os.path.exists('wannier90.win'):
        raise FileNotFoundError("wannier90.win not found in current directory")
    if not os.path.exists('wannier90.eig'):
        raise FileNotFoundError("wannier90.eig not found in current directory")
    if not os.path.exists('wannier90_hr.dat'):
        raise FileNotFoundError("wannier90_hr.dat not found in current directory")
    if not os.path.exists('FERMI_ENERGY.in'):
        raise FileNotFoundError("FERMI_ENERGY.in not found in current directory")
    
    # 读取数据并生成文件
    lattice, num_atoms,atoms = read_poscar('POSCAR')
    projections = read_wannier_win('wannier90.win')
    fermi,num_ocps = get_fermi_occupation("wannier90.eig","wannier90.win","FERMI_ENERGY.in")
    soc_val=1 if soc else 0
    write_wt_in(lattice, num_atoms,atoms, projections,fermi,num_ocps,kpath,soc=soc_val)
    write_job(main_path)
    print("Successfully generated wt.in file!")

if __name__ == "__main__":
    main()
