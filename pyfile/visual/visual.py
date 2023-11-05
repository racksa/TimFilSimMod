import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import myIo
import util
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import matplotlib.animation as animation
import configparser
import math
import sys
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

class VISUAL:

    def __init__(self):
        self.globals_name = 'globals.ini'
        self.dir = "/home/clustor2/ma/h/hs2216/20231027/"
        # self.dir = "data/expr_sims/20231027/"
        self.dir = "data/expr_sims/20231105/"
        self.pars_list = {
                     "nswim": [],
                     "nseg": [],
                     "nfil": [],
                     "nblob": [],
                     "ar": [],
                     "spring_factor": []}
        self.video = False
        self.interpolate = False
        self.angle = False
        self.output_to_fcm = False
        self.output_to_superpunto = True
        self.periodic = False
        self.big_sphere = True
        self.check_overlap = False

        self.plot_end_frame_setting = 30000
        self.frames = 90

        self.plot_end_frame = self.plot_end_frame_setting

        self.select_elst = 0

        self.Lx = 1000
        self.Ly = 1000
        self.Lz = 1000
        
        self.index = 0

    def write_data(self, x, r, filename, box=True, center=True, superpunto=True, color=0):
        plot_x = x.copy()
        if(box):
            plot_x[0] = util.box(plot_x[0], self.Lx)
            plot_x[1] = util.box(plot_x[1], self.Ly)
            plot_x[2] = util.box(plot_x[2], self.Lz)
            if(center):
                plot_x[0] -= 0.5*self.Lx
                plot_x[1] -= 0.5*self.Ly
        if superpunto:
            plot_x[1] = -plot_x[1]
            plot_x[2] = -plot_x[2]
            myIo.write_line(str(plot_x[0]) + ' ' +\
                            str(plot_x[1]) + ' ' +\
                            str(plot_x[2]) + ' ' +\
                            str(r) + ' ' +\
                            str(color),
                            filename)
        else:
            myIo.write_line(str(plot_x[0]) + ' ' +\
                            str(plot_x[1]) + ' ' +\
                            str(plot_x[2]),
                            filename)
        

    def read_rules(self):
        sim = configparser.ConfigParser()
        try:
            sim.read(self.dir+"rules.ini")
            num_fil = len(np.unique([float(s) for s in sim["Parameter list"]['nfil'].split(', ')]))
            num_ar = len(np.unique([float(s) for s in sim["Parameter list"]['ar'].split(', ')]))
            num_elst = len(np.unique([float(s) for s in sim["Parameter list"]['spring_factor'].split(', ')]))
            num_per_elst = int(num_fil*num_ar)
            select_elst = min(num_elst-1, self.select_elst)

            for key, value in self.pars_list.items():
                if(key in sim["Parameter list"]):
                    # self.pars_list[key] = [float(x) for x in sim["Parameter list"][key].split(', ')][0::num_elst]
                    self.pars_list[key] = [float(x) for x in sim["Parameter list"][key].split(', ')][num_per_elst*select_elst:num_per_elst*(select_elst+1)]
            self.num_sim = len(self.pars_list["nfil"])
        except:
            print("WARNING: " + self.dir + "rules.ini not found.")

    def select_sim(self):
        if(self.index>len(self.pars_list['nfil'])):
            self.index = len(self.pars_list['nfil'])-1
            print(f'Index out of range. Using the last sim: {self.index}')

        self.nseg = int(self.pars_list['nseg'][self.index])
        self.nswim = 1
        self.nfil = int(self.pars_list['nfil'][self.index])
        self.nblob = int(self.pars_list['nblob'][self.index])
        self.ar = self.pars_list['ar'][self.index]
        self.spring_factor = self.pars_list['spring_factor'][self.index]
        self.N = int(self.nswim*(self.nfil*self.nseg + self.nblob))

        self.simName = self.dir + f"ciliate_{self.nfil:.0f}fil_{self.nblob:.0f}blob_{self.ar:.2f}R_{self.spring_factor:.2f}torsion"
        self.fil_references = myIo.read_fil_references(self.simName + '_fil_references.dat')

        self.pars = myIo.read_pars(self.simName + '.par')
        self.radius = 0.5*self.ar*self.pars['FIL_LENGTH']
        self.dt = self.pars['DT']*self.pars['PLOT_FREQUENCY_IN_STEPS']
        if(not 'PRESCRIBED_CILIA' in self.pars):
            self.pars['PRESCRIBED_CILIA'] = 0
        if(self.pars['NBLOB']>0):
            self.blob_references = myIo.read_blob_references(self.simName + '_blob_references.dat')
        if(self.pars['NFIL']>0):
            self.fil_references = myIo.read_fil_references(self.simName + '_fil_references.dat')

        self.plot_end_frame = min(self.plot_end_frame_setting, sum(1 for line in open(self.simName + '_body_states.dat')))
        self.plot_start_frame = max(0, self.plot_end_frame-self.frames)
        self.frames = self.plot_end_frame - self.plot_start_frame
        self.plot_interval = 1
        print(f'index={self.index} file={self.simName}')
        
    def plot(self):
        self.select_sim()
        ax = plt.figure().add_subplot(projection='3d')

        superpuntoDatafileName = self.simName + "_superpunto.dat"
        myIo.clean_file(superpuntoDatafileName)

        seg_states_f = open(self.simName + '_seg_states.dat', "r")
        body_states_f = open(self.simName + '_body_states.dat', "r")
        # if(self.pars['NBLOB']>0):
        #     blob_forces_f = open(self.simName + '_blob_forces.dat', "r")
        # seg_forces_f = open(self.simName + '_seg_forces.dat', "r")
        if (self.pars['PRESCRIBED_CILIA'] == 1):
            fil_phases_f = open(self.simName + '_filament_phases.dat', "r")
            # fil_angles_f = open(self.simName + '_filament_shape_rotation_angles.dat', "r")

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            if(self.check_overlap):
                segs_list = np.zeros((int(self.nfil*self.nseg), 3))
                blobs_list = np.zeros((int(self.nblob), 3))
                

            body_states_str = body_states_f.readline()
            if(self.pars['NFIL']>0):
                seg_states_str = seg_states_f.readline()
                if (self.pars['PRESCRIBED_CILIA'] == 1):
                    fil_phases_str = fil_phases_f.readline()
                    # fil_angles_str = fil_angles_f.readline()

                    fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                    fil_phases = util.box(fil_phases, 2*np.pi)

            if(i%self.plot_interval==0 and i>=self.plot_start_frame):
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                if(self.pars['NFIL']>0):
                    seg_states = np.array(seg_states_str.split()[1:], dtype=float)                
                
                myIo.write_line('#', superpuntoDatafileName)
                for swim in range(int(self.pars['NSWIM'])):
                    body_pos = body_states[7*swim : 7*swim+3]
                    R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                    R = np.linalg.inv(R)
                    R = np.eye(3)
                    if(not self.big_sphere or self.check_overlap):
                        for blob in range(int(self.pars['NBLOB'])):
                            blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                            if(self.check_overlap):
                                blobs_list[blob] = blob_x, blob_y, blob_z
                            elif(not self.big_sphere):
                                self.write_data([blob_x, blob_y, blob_z], float(self.pars['RBLOB']), superpuntoDatafileName, self.periodic, color=16777215)
                            
                    if(self.big_sphere):
                        self.write_data(body_pos, self.radius, superpuntoDatafileName, self.periodic, color=16777215)

                    for fil in range(int(self.pars['NFIL'])):
                        fil_color = int("000000", base=16)
                        # Robot arm to find segment position (Ignored plane rotation!)
                        if (self.pars['PRESCRIBED_CILIA'] == 0):
                            fil_i = int(4*fil*self.pars['NSEG'])
                            fil_base_x, fil_base_y, fil_base_z = body_pos + np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                            old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z])
                        elif (self.pars['PRESCRIBED_CILIA'] == 1):
                            fil_i = int(3*fil*self.pars['NSEG'])
                            old_seg_pos = seg_states[fil_i : fil_i+3]
                            
                            # WRITE A FUNCTION FOR THIS!!
                            cmap_name = 'hsv'
                            cmap = plt.get_cmap(cmap_name)
                            rgb_color = cmap(fil_phases[fil]/(2*np.pi))[:3]  # Get the RGB color tuple
                            rgb_hex = mcolors.rgb2hex(rgb_color)[1:]  # Convert RGB to BGR hexadecimal format
                            bgr_hex = rgb_hex[4:]+rgb_hex[2:4]+rgb_hex[:2]
                            fil_color = int(bgr_hex, base=16)
                            # print("\n", bgr_hex, fil_color, "\t")
                        if(self.check_overlap):
                            segs_list[fil*self.nseg] = old_seg_pos
                        else:
                            self.write_data(old_seg_pos, float(self.pars['RSEG']), superpuntoDatafileName, self.periodic, True, True, color=fil_color)

                        for seg in range(1, int(self.pars['NSEG'])):
                            if (self.pars['PRESCRIBED_CILIA'] == 0):
                                q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                                q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                                
                                t1 = util.find_t(q1)
                                t2 = util.find_t(q2)
                                
                                seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                                old_seg_pos = seg_pos
                            elif (self.pars['PRESCRIBED_CILIA'] == 1):
                                seg_pos = seg_states[fil_i+3*(seg) : fil_i+3*(seg+1)] 
                            if(self.check_overlap):
                                segs_list[fil*self.nseg + seg] = seg_pos
                            else:
                                self.write_data(seg_pos, float(self.pars['RSEG']), superpuntoDatafileName, self.periodic, True, True, color=fil_color)
                            
                if(self.check_overlap):
                    threshold = 1.00
                    colliding_indices, colliding_particles = util.label_colliding_particles_with_3d_cell_list(segs_list, 5, threshold*float(self.pars['RSEG']))
                    
                    self.write_data(body_pos, self.radius, superpuntoDatafileName, self.periodic, color=16777215)
                    print(f'Overlapping case at threshold {threshold} = {len(colliding_indices)}')

                    for i, pos in enumerate(segs_list):
                        fil_color = 16777215
                        if(i in colliding_indices):
                            fil_color = 255
                        self.write_data(pos, float(self.pars['RSEG']), superpuntoDatafileName, self.periodic, True, True, color=fil_color)
                        

## Filaments
    def plot_fil(self):
        self.select_sim()

        seg_states_f = open(self.simName + '_seg_states.dat', "r")
        body_states_f = open(self.simName + '_body_states.dat', "r")
        
        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_proj_type('ortho')
        # ax.set_proj_type('persp', 0.05)  # FOV = 157.4 deg
        # ax.view_init(elev=5., azim=45)
        # ax.dist=20
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        # ax.axis('off')
        # ax.grid(False)


        def animation_func(t):
            print(t)
            ax.cla()

            ax.set_xlim(-100, 300)
            ax.set_ylim(-100, 300)
            ax.set_zlim(-100, 300)

            seg_states_str = seg_states_f.readline()

            seg_states = np.array(seg_states_str.split()[1:], dtype=float)

            # Robot arm to find segment position (Ignored plane rotation!)
            for fil in range(self.nfil):
                fil_data = np.zeros((self.nseg, 3))
                fil_i = int(3*fil*self.nseg)
                fil_data[0] = seg_states[fil_i : fil_i+3]

                for seg in range(1, self.nseg):
                    seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                    fil_data[seg] = seg_pos
                ax.plot(fil_data[:,0], fil_data[:,1], fil_data[:,2], c='black', zorder = 100)

        if(self.video):
            plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
            ani = animation.FuncAnimation(fig, animation_func, frames=500, interval=10, repeat=False)
            plt.show()
            # FFwriter = animation.FFMpegWriter(fps=10)
            # ani.save(f'fig/ciliate_{nfil}fil_anim.mp4', writer=FFwriter)
        else:
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                if(i==self.plot_end_frame-1):
                    animation_func(i)
                else:
                    seg_states_str = seg_states_f.readline()
                
            plt.savefig(f'fig/ciliate_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
            plt.show()


## Ciliates
# Single sim
    def phase(self):
        self.select_sim()
        
        fil_phases_f = open(self.simName + '_filament_phases.dat', "r")
        fil_angles_f = open(self.simName + '_filament_shape_rotation_angles.dat', "r")

        # Plotting
        colormap = 'cividis'
        colormap = 'twilight_shifted'

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fil_references_sphpolar = np.zeros((self.nfil,3))

        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        norm = Normalize(vmin=0, vmax=2*np.pi)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.ax.set_yticks(np.linspace(0, 2*np.pi, 7), ['0', 'π/3', '2π/3', 'π', '4π/3', '5π/3', '2π'])
        cbar.set_label(r"phase")

        ax.set_ylabel(r"$\theta$")
        ax.set_xlabel(r"$\phi$")
        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(0, np.pi)
        ax.set_xticks(np.linspace(-np.pi, np.pi, 5), ['-π', '-π/2', '0', 'π/2', 'π'])
        ax.set_yticks(np.linspace(0, np.pi, 5), ['0', 'π/4', 'π/2', '3π/4', 'π'])

        global frame
        frame = 0
        import scipy.interpolate

        def animation_func(t):
            global frame

            ax.cla()
            fil_phases_str = fil_phases_f.readline()
            fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
            fil_phases = util.box(fil_phases, 2*np.pi)

            fil_angles_str = fil_angles_f.readline()
            fil_angles = np.array(fil_angles_str.split()[1:], dtype=float)

            for i in range(self.nfil):
                fil_references_sphpolar[i] = util.cartesian_to_spherical(self.fil_references[3*i: 3*i+3])

            if self.angle:
                variables = fil_angles
            else:
                variables = fil_phases

            # Interpolation
            if (self.interpolate):
                n1, n2 = 100, 100
                azim_grid = np.linspace(-np.pi, np.pi, n1)
                polar_grid = np.linspace(0, np.pi, n2)
                xx, yy = np.meshgrid(azim_grid, polar_grid)
                zz = scipy.interpolate.griddata((fil_references_sphpolar[:,1],fil_references_sphpolar[:,2]), variables, (xx, yy), method='cubic')
                ax.scatter(xx, yy, c=zz, cmap=colormap)
            else:
            # Individual filaments
                ax.scatter(fil_references_sphpolar[:,1], fil_references_sphpolar[:,2], c=variables, cmap=colormap)

            frame += 1

        if(self.video):
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                if(i>=self.plot_start_frame):
                    frame = i
                    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
                    ani = animation.FuncAnimation(fig, animation_func, frames=500, interval=10, repeat=False)
                    plt.show()
                    FFwriter = animation.FFMpegWriter(fps=10)
                    ani.save(f'fig/fil_phase_{self.nfil}fil_anim.mp4', writer=FFwriter)
                    break
                else:
                    fil_phases_str = fil_phases_f.readline()
        else:
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                if(i==self.plot_end_frame-1):
                    animation_func(i)
                else:
                    fil_phases_str = fil_phases_f.readline()

            plt.savefig(f'fig/fil_phase_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
            plt.show()

    def order_parameter(self):
        self.select_sim()

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1,1,1)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(1,1,1)
        fil_references_sphpolar = np.zeros((self.nfil,3))
        
        fil_phases_f = open(self.simName + '_filament_phases.dat', "r")
        fil_angles_f = open(self.simName + '_filament_shape_rotation_angles.dat', "r")


        nfil = self.nfil
        n_snapshots = min(30, self.plot_end_frame)
        start = self.plot_end_frame - n_snapshots
        X = np.zeros((nfil, n_snapshots))
        time_array = np.arange(n_snapshots)

        fil_references_sphpolar = np.zeros((nfil,3))
        for fil in range(nfil):
            fil_references_sphpolar[fil] = util.cartesian_to_spherical(self.fil_references[3*fil: 3*fil+3])
        azim_array = fil_references_sphpolar[:,1]
        polar_array = fil_references_sphpolar[:,2]
        sorted_indices = np.argsort(azim_array)
        azim_array_sorted = azim_array[sorted_indices]
        polar_array_sorted = polar_array[sorted_indices]
        

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            fil_phases_str = fil_phases_f.readline()
            fil_angles_str = fil_angles_f.readline()
            
            if(i>=start):
                fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                fil_phases_sorted = fil_phases[sorted_indices]
                fil_phases_sorted = util.box(fil_phases_sorted, 2*np.pi)

                fil_angles = np.array(fil_angles_str.split()[1:], dtype=float)
                fil_angles_sorted = fil_angles[sorted_indices]

                # fil_phases_sorted = np.sin(fil_phases_sorted)

                X[:,i-start] = fil_phases_sorted[:nfil]
                # X[:,i-start] = fil_angles_sorted[:nfil]
        
        
        
        phi_dot = X[:,1:] - X[:,:-1]
        for i in range(5):
            ax.plot(time_array[1:], phi_dot[i,:])
        print(max(phi_dot[:,-1]), min(phi_dot[:,-1]))
        ax.set_xlabel('time')

        # phi_dot = phi_dot[:,0]

        # phi_dot_mean = np.mean()

        r = np.abs(np.sum(np.exp(phi_dot*1j), 0)/self.nfil)
        print(f"index={self.index} r={np.mean(r)}")
        # print(np.exp(phi_dot*1j))
        # print(np.sum(np.exp(phi_dot*1j), 0)/self.nfil)

        # print(self.nfil, np.shape(r))
        
        # print(r)
        ax2.scatter(azim_array_sorted, polar_array_sorted, c = phi_dot[:,-1])

        # ax3.plot(time_array[1:], r)

        fig.savefig(f'fig/fil_order_parameter_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()
        
    def eckert(self):
        R = 1
        phi0 = np.pi
        from scipy import optimize
        def find_k(k, theta):
            return k + np.sin(k)*np.cos(k) + 2*np.sin(k) - (2 + np.pi/2)*np.sin(theta)

        def eckert_projection(theta, phi):
            sign = 1
            if theta >= np.pi/2:
                theta = (np.pi/2 - (theta - np.pi/2))
                sign = -1
            
            k=optimize.newton(find_k, 1, args=(theta,))

            x = 0.4222382*R*(phi-phi0)*(1+np.cos(k))
            y = 1.3265004*R*np.sin(k)*sign
            return x, y
    
        self.select_sim()
        
        fil_phases_f = open(self.simName + '_filament_phases.dat', "r")
        fil_angles_f = open(self.simName + '_filament_shape_rotation_angles.dat', "r")

        # Plotting
        colormap = 'cividis'
        colormap = 'twilight_shifted'

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fil_references_sphpolar = np.zeros((self.nfil,3))

        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        norm = Normalize(vmin=0, vmax=2*np.pi)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.ax.set_yticks(np.linspace(0, 2*np.pi, 7), ['0', 'π/3', '2π/3', 'π', '4π/3', '5π/3', '2π'])
        cbar.set_label(r"phase")

        ax.set_ylabel(r"$\theta$")
        ax.set_xlabel(r"$\phi$")
        # ax.set_xlim(-np.pi, np.pi)
        # ax.set_ylim(0, np.pi)
        ax.set_xticks(np.linspace(-np.pi, np.pi, 5), ['-π', '-π/2', '0', 'π/2', 'π'])
        ax.set_yticks(np.linspace(0, np.pi, 5), ['0', 'π/4', 'π/2', '3π/4', 'π'])

        global frame
        frame = 0

        def animation_func(t):
            global frame

            ax.cla()
            fil_phases_str = fil_phases_f.readline()

            fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
            fil_phases = util.box(fil_phases, 2*np.pi)
            for i in range(self.nfil):
                fil_references_sphpolar[i] = util.cartesian_to_spherical(self.fil_references[3*i: 3*i+3])
                
            projected_points = [eckert_projection(theta, phi) for theta, phi in zip(fil_references_sphpolar[:,2], fil_references_sphpolar[:,1])]
            projected_x, projected_y = zip(*projected_points)
            
            ax.scatter(projected_x, projected_y, c=fil_phases, cmap=colormap)
            frame += 1

        if(self.video):
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                if(i>=self.plot_start_frame):
                    frame = i
                    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
                    ani = animation.FuncAnimation(fig, animation_func, frames=500, interval=10, repeat=False)
                    plt.show()
                    FFwriter = animation.FFMpegWriter(fps=10)
                    ani.save(f'fig/fil_phase_{self.nfil}fil_anim.mp4', writer=FFwriter)
                    break
                else:
                    fil_phases_str = fil_phases_f.readline()
        else:
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                if(i==self.plot_end_frame-1):
                    animation_func(i)
                else:
                    fil_phases_str = fil_phases_f.readline()

            plt.savefig(f'fig/fil_phase_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
            plt.show()
        

    def ciliate(self):
        self.select_sim()

        seg_states_f = open(self.simName + '_seg_states.dat', "r")
        body_states_f = open(self.simName + '_body_states.dat', "r")

        # Create the sphere data points
        num_points = 300
        u = np.linspace(0, 2 * np.pi, num_points)
        v = np.linspace(0, np.pi, num_points)
        x = self.radius * np.outer(np.cos(u), np.sin(v))
        y = self.radius * np.outer(np.sin(u), np.sin(v))
        z = self.radius * np.outer(np.ones(np.size(u)), np.cos(v))

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_proj_type('ortho')
        # ax.set_proj_type('persp', 0.05)  # FOV = 157.4 deg
        # ax.view_init(elev=5., azim=45)
        # ax.dist=20
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        # ax.axis('off')
        # ax.grid(False)

        def animation_func(t):
            print(t)
            ax.cla()

            # ax.set_xlim(-100, 100)
            # ax.set_ylim(-100, 100)
            # ax.set_zlim(-100, 100)

            body_states_str = body_states_f.readline()
            seg_states_str = seg_states_f.readline()

            body_states = np.array(body_states_str.split()[1:], dtype=float)
            seg_states = np.array(seg_states_str.split()[1:], dtype=float)
            
            for swim in range(self.nswim):
                # blob_data = np.zeros((int(self.pars['NBLOB']), 3))
                body_pos = body_states[7*swim : 7*swim+3]
                # R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                # for blob in range(int(self.pars['NBLOB'])):
                #     blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                #     ax.scatter(blob_x, blob_y, blob_z)

                # Plot the sphere
                ax.plot_surface(x+body_pos[0], y+body_pos[1], z+body_pos[2], color='grey', alpha=0.5)

                # Robot arm to find segment position (Ignored plane rotation!)
                for fil in range(self.nfil):
                    fil_data = np.zeros((self.nseg, 3))
                    fil_i = int(3*fil*self.nseg)
                    fil_data[0] = seg_states[fil_i : fil_i+3]

                    for seg in range(1, self.nseg):
                        seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                        fil_data[seg] = seg_pos
                    ax.plot(fil_data[:,0], fil_data[:,1], fil_data[:,2], c='black', zorder = 100)

        if(self.video):
            plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
            ani = animation.FuncAnimation(fig, animation_func, frames=500, interval=10, repeat=False)
            plt.show()
            # FFwriter = animation.FFMpegWriter(fps=10)
            # ani.save(f'fig/ciliate_{nfil}fil_anim.mp4', writer=FFwriter)
        else:
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                if(i==self.plot_end_frame-1):
                    animation_func(i)
                else:
                    body_states_str = body_states_f.readline()
                    seg_states_str = seg_states_f.readline()
                
            plt.savefig(f'fig/ciliate_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
            plt.show()

    def ciliate_traj(self):
        self.select_sim()

        body_states_f = open(self.simName + '_body_states.dat', "r")
        time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
        body_pos_array = np.zeros((len(time_array), 3))

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            body_states_str = body_states_f.readline()

            if(i>=self.plot_start_frame):
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                body_pos_array[i-self.plot_start_frame] = body_states[0 : 3]

        ax.plot(body_pos_array[:,0], body_pos_array[:,1], body_pos_array[:,2])
        plt.savefig(f'fig/ciliate_traj_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def ciliate_speed(self):
        self.select_sim()

        # seg_states_f = open(self.simName + '_seg_states.dat', "r")
        body_states_f = open(self.simName + '_body_states.dat', "r")
        body_vels_f = open(self.simName + '_body_vels.dat', "r")

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
        
        body_pos_array = np.zeros((len(time_array), 3))
        body_vel_array = np.zeros((len(time_array), 6))
        body_speed_array = np.zeros(len(time_array))

        # pos = np.zeros(3)

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            body_states_str = body_states_f.readline()
            body_vels_str = body_vels_f.readline()

            if(i>=self.plot_start_frame):

                body_states = np.array(body_states_str.split()[1:], dtype=float)
                body_vels = np.array(body_vels_str.split(), dtype=float)

                body_pos_array[i-self.plot_start_frame] = body_states[0:3]
                body_vel_array[i-self.plot_start_frame] = body_vels
                body_speed_array[i-self.plot_start_frame] = np.sqrt(np.sum(body_vels[0:3]*body_vels[0:3], 0))

            # pos += body_vels[0:3]
            # body_vel_array[i][0:3] = pos*self.dt

        # for i in range(len(body_vel_array[0])):
        #     # ax.plot(time_array, body_pos_array[:,i])
        #     ax.plot(time_array, body_vel_array[:,i])
        ax.plot(time_array, body_speed_array)
        plt.savefig(f'fig/ciliate_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def ciliate_forcing(self):
        self.select_sim()

        seg_forces_f = open(self.simName + '_seg_forces.dat', "r")
        blob_forces_f = open(self.simName + '_blob_forces.dat', "r")
        

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
        body_force_array = np.zeros((len(time_array), 3))

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            seg_forces_str = seg_forces_f.readline()
            blob_forces_str = blob_forces_f.readline()

            if(i>=self.plot_start_frame):

                seg_forces = np.array(seg_forces_str.split()[1:], dtype=float)
                blob_forces= np.array(blob_forces_str.split()[1:], dtype=float)

                seg_forces = np.reshape(seg_forces, (int(self.pars['NSEG']*self.pars['NFIL']), 6))
                blob_forces = np.reshape(blob_forces, (int(self.pars['NBLOB']), 3))

                body_force_array[i-self.plot_start_frame] = np.sum(blob_forces, axis=0) + np.sum(seg_forces[:,0:3], axis=0) 

        labels=[r'$\lambda_x$', r'$\lambda_y$', r'$\lambda_z$']
        for i in range(len(body_force_array[0])):
            ax.plot(time_array, body_force_array[:,i], label=labels[i])
        ax.set_xlabel('Time step')
        ax.set_ylabel('Force')
        plt.legend()
        plt.savefig(f'fig/ciliate_forcing_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def ciliate_dissipation(self):
        self.select_sim()

        seg_forces_f = open(self.simName + '_seg_forces.dat', "r")
        seg_vels_f = open(self.simName + '_seg_vels.dat', "r")
        blob_forces_f = open(self.simName + '_blob_forces.dat', "r")
        blob_references_f = open(self.simName + '_blob_references.dat', "r")
        body_vels_f = open(self.simName + '_body_vels.dat', "r")

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
        dissipation_array = np.zeros(self.frames)

        blob_references_str = blob_references_f.readline()
        blob_references= np.array(blob_references_str.split(), dtype=float)
        blob_references = np.reshape(blob_references, (int(self.pars['NBLOB']), 3))

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            seg_forces_str = seg_forces_f.readline()
            seg_vels_str = seg_vels_f.readline()
            blob_forces_str = blob_forces_f.readline()
            body_vels_str = body_vels_f.readline()

            if(i>=self.plot_start_frame):

                seg_forces = np.array(seg_forces_str.split()[1:], dtype=float)
                seg_vels = np.array(seg_vels_str.split()[1:], dtype=float)
                blob_forces= np.array(blob_forces_str.split()[1:], dtype=float)
                body_vels= np.array(body_vels_str.split(), dtype=float)

                seg_forces = np.reshape(seg_forces, (int(self.pars['NSEG']*self.pars['NFIL']), 6))
                seg_vels = np.reshape(seg_vels, (int(self.pars['NSEG']*self.pars['NFIL']), 6))
                blob_forces = np.reshape(blob_forces, (int(self.pars['NBLOB']), 3))
                body_vels_tile = np.tile(body_vels, (int(self.pars['NBLOB']), 1))
                blob_vels = body_vels_tile[:, 0:3] + np.cross(body_vels_tile[:, 3:6], blob_references)

                dissipation_array[i-self.plot_start_frame] = np.sum(blob_forces * blob_vels) + np.sum(seg_forces * seg_vels)

        ax.plot(time_array, dissipation_array)
        ax.set_xlabel('Time step')
        ax.set_ylabel('Dissipation')
        plt.savefig(f'fig/ciliate_dissipation_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def ciliate_dmd(self):
        self.select_sim()
        
        fil_phases_f = open(self.simName + '_filament_phases.dat', "r")


        nfil = self.nfil
        n_snapshots = min(30, self.plot_end_frame)
        start = self.plot_end_frame - n_snapshots
        X = np.zeros((nfil, n_snapshots))
        r = min(nfil, n_snapshots-1)
        print(f'r={r}')
        dt = 1./30
        coeffs = np.zeros((r, n_snapshots), dtype=complex)

        fil_references_sphpolar = np.zeros((nfil,3))
        for fil in range(nfil):
            fil_references_sphpolar[fil] = util.cartesian_to_spherical(self.fil_references[3*fil: 3*fil+3])
        azim_array = fil_references_sphpolar[:,1]
        sorted_indices = np.argsort(azim_array)
        azim_array_sorted = azim_array[sorted_indices]

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            fil_phases_str = fil_phases_f.readline()
            
            if(i>=start):
                fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                fil_phases_sorted = fil_phases[sorted_indices]
                fil_phases_sorted = util.box(fil_phases_sorted, 2*np.pi)

                fil_phases_sorted = np.sin(fil_phases_sorted)

                X[:,i-start] = fil_phases_sorted[:nfil]
                
        X1 = X[:, :-1]
        X2 = X[:, 1:]

        U, sigma, V = np.linalg.svd(X1, full_matrices=False)
        Sigma = np.diag(sigma)
        V = V.conj().T
        U = U[:, :r]
        Sigma = Sigma[:r, :r]
        V = V[:, :r]
        


        A_tilde = U.conj().T @ X2 @ V @ np.linalg.inv(Sigma)
        D, W = np.linalg.eig(A_tilde)
        omega = np.log(D)/dt
        phi = X2 @ V @ np.linalg.inv(Sigma) @ W
        b = np.linalg.pinv(phi) @ X1[:,0]

        for t in range(n_snapshots):
            coeffs[:, t] = np.exp(omega * t*dt) * b

        X_dmd = (phi @ coeffs).real

            
        # print(np.shape(D), np.shape(W), np.shape(A_tilde))
        # print(np.shape(b))
        # print(np.shape(omega))
        # print(phi)
        # print(np.shape(np.exp(omega * dt)))

        inspected_snapshots = np.array([0, 1, 2])
        modes = np.abs(b).argsort()[-4:][::-1]
        # modes = [0,1,2,3]
    
        # Plotting
        fig, axs = plt.subplots(len(modes), sharex=True, sharey=True)
        axs_flat = axs.ravel()
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1,1,1)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(1,1,1)
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(1,1,1)

        print(np.shape(phi))
        
        for ind, mode in enumerate(modes):
            axs[ind].plot(azim_array_sorted, phi[:, mode].real, label=f'real')
            axs[ind].plot(azim_array_sorted, phi[:, mode].imag, label=f'imag')

            # axs[ind].plot(coeffs[:, mode].real, label=f'real')
            # axs[ind].plot(coeffs[:, mode].imag, label=f'imag')

            # axs[ind].set_xlabel('Azimuth angle')
            # axs[ind].set_ylabel(r'$\phi$')
            axs[ind].set_title(f'mode={mode}')
            axs[ind].legend()
        axs[-1].set_xlabel("Azimuthal position")

        ax2.plot(np.abs(b))
        # for time in range(n_snapshots):
        #     ax2.plot(np.abs(coeffs[:, time]))
        ax2.set_xlabel(r'Mode')
        ax2.set_ylabel(r'Magnitude(b)')

        for i in inspected_snapshots:
            ax3.plot(X[:, i], c='b', marker='+')
            ax3.plot(X_dmd[:, i], c='r', marker='x')
        ax3.set_xlabel("fil")
        ax3.set_ylabel(r"sin(phase)")

        ax4.imshow(X)
        # ax4.plot(sigma)
        ax4.set_xlabel("t")
        ax4.set_ylabel(r"fil i")
        

        fig.savefig(f'fig/fil_dmd_modes_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        # fig2.savefig(f'fig/fil_svd_cumsum_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        # fig3.savefig(f'fig/fil_svd_modes_{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def ciliate_svd(self):
        self.select_sim()
        
        fil_phases_f = open(self.simName + '_filament_phases.dat', "r")

        nfil = self.nfil
        n_snapshots = min(30, self.plot_end_frame)
        start = self.plot_end_frame - n_snapshots
        X = np.zeros((nfil, n_snapshots))

        fil_references_sphpolar = np.zeros((nfil,3))
        for fil in range(nfil):
            fil_references_sphpolar[fil] = util.cartesian_to_spherical(self.fil_references[3*fil: 3*fil+3])
        azim_array = fil_references_sphpolar[:,1]
        polar_array = fil_references_sphpolar[:,2]
        sorted_indices = np.argsort(azim_array)
        azim_array_sorted = azim_array[sorted_indices]
        polar_array_sorted = polar_array[sorted_indices]
        

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
            fil_phases_str = fil_phases_f.readline()
            
            if(i>=start):
                fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                fil_phases_sorted = fil_phases[sorted_indices]
                fil_phases_sorted = util.box(fil_phases_sorted, 2*np.pi)

                fil_phases_sorted = np.sin(fil_phases_sorted)

                X[:,i-start] = fil_phases_sorted[:nfil]

        U, sigma, V = np.linalg.svd(X, full_matrices=False)
        Sigma = np.diag(sigma)
        # V = V.conj().T
        # U = U[:, :r]
        # Sigma = Sigma[:r, :r]
        # V = V[:, :r]

        # res = np.linalg.svd(X)
        # svd_diag = np.zeros(np.shape(X))
        # diag = np.diag(res[1])
        # svd_diag[:diag.shape[0], :diag.shape[1]] = diag

        # print(np.shape(U),np.shape(Sigma), np.shape(V))

        pc = U @ Sigma
        pa = V        
        
        num_fil = 4
        num_mode = 2
        reduced = pc @ pa

        

        # Plotting
        colormap = 'twilight_shifted'
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        norm = Normalize(vmin=0, vmax=2*np.pi)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1,1,1)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(1,1,1)
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(1,1,1)
        fig5 = plt.figure()
        ax5 = fig5.add_subplot(1,1,1)
        fig6 = plt.figure()
        ax6 = fig6.add_subplot(1,1,1)
        fig7 = plt.figure()
        ax7 = fig7.add_subplot(1,1,1)
        fig8 = plt.figure()
        ax8 = fig8.add_subplot(1,1,1)

        # Signal X
        ax.imshow(X)
        ax.set_xlabel('t')
        ax.set_ylabel('Fil index (sorted by azimuth angle)')

        # Signal of the first snapshot
        ax2.scatter(azim_array_sorted, X[:,0])
        ax2.set_xlabel('Azimuth position')
        ax2.set_ylabel('Phase, t=0')

        # Eigenvalue
        ax3.scatter(np.arange(0, len(sigma), 1), sigma, marker='*')
        ax3.set_xlabel(r'Mode')
        ax3.set_ylabel(r'$Eigenvalue$')

        # POD spatial modes
        for i in range(2):
            ax4.plot(U[:,i], label=f'Mode {i}')
            ax4.set_xlabel('x')
            ax4.set_ylabel('f(x, t=0)')
            ax4.legend()

        # POD temporal modes
        for i in range(4):
            ax5.plot(V[i, :], label=f'Mode {i}')
            ax5.set_xlabel('t')
            ax5.set_ylabel('f(x, t=0)')
            ax5.legend()

        # Phase of the first snapshot
        ax6.scatter(azim_array_sorted, polar_array_sorted, c=X[:,0], cmap=colormap)
        ax6.set_xlabel(r"Azimuth position")
        ax6.set_ylabel(r"Polar position")

        # Interpolated phase of the first snapshot
        n1, n2 = 100, 100
        azim_grid = np.linspace(-np.pi, np.pi, n1)
        polar_grid = np.linspace(0, np.pi, n2)
        xx, yy = np.meshgrid(azim_grid, polar_grid)
        import scipy.interpolate
        zz = scipy.interpolate.griddata((azim_array_sorted, polar_array_sorted), X[:,20], (xx, yy), method='cubic')
        ax7.scatter(xx, yy, c=zz, cmap=colormap)
        ax7.set_xlabel(r"Azimuth position")
        ax7.set_ylabel(r"Polar position")
        
        # Interpolated phase of reconstruction using first mode
        uu = scipy.interpolate.griddata((azim_array_sorted, polar_array_sorted), U[:,1], (xx, yy), method='cubic')
        ax8.scatter(xx, yy, c=uu, cmap=colormap)
        ax8.set_xlabel(r"Azimuth position")
        ax8.set_ylabel(r"Polar position")



        # ax.scatter(azim_array_sorted, fil_phases_sorted)
        # for fil in range(num_fil):
        #     abs_pc = np.abs(pc[fil][:nfil])
        #     ax2.plot(np.cumsum(abs_pc)/np.sum(abs_pc), label=f'fil {fil}')
        # ax.set_xlabel('Azimuth angle')
        # ax.set_ylabel(r'$\phi$')

        # ax2.set_xlabel('Mode')
        # ax2.set_ylabel('Accumulated |weight| fraction')
        # ax2.legend()

        
        
        # for i in range(6):
        #     ax4.plot(U[i], c='b')
        #     ax4.plot(reduced[i], c='r')

        # ax4.plot(0, c='b', label='Original' )
        # ax4.plot(0, c='r', label=f'Using {num_mode} modes')
        # ax4.set_xlabel('Time step')
        # ax4.set_ylabel(r'$\phi$')
        # ax4.legend()
        

        fig.savefig(f'fig/fil_svd_signal_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        fig2.savefig(f'fig/fil_svd_initial_snapshot_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        fig3.savefig(f'fig/fil_svd_eigenvalues_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        fig4.savefig(f'fig/fil_svd_spatial_modes_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        fig5.savefig(f'fig/fil_svd_temporal_modes_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        fig6.savefig(f'fig/fil_svd_phase_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        fig7.savefig(f'fig/fil_svd_interpolated_phase_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        fig8.savefig(f'fig/fil_svd_first_mode_reconstruction_{self.nfil}fil_{self.index}.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def timing(self):
        self.select_sim()

        with open(self.simName + '_time.dat', 'r') as file:
            plot_time_frame = len(file.readlines())
        timings_f = open(self.simName + '_time.dat', "r")

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
        time_array = np.arange(0, plot_time_frame)
        
        time_end_frame = plot_time_frame
        time_start_frame = max(0, time_end_frame-100)
        timings_array = np.zeros((plot_time_frame, 7))
        total_timing_array = np.zeros(plot_time_frame)

        for i in range(plot_time_frame):
            print(" frame ", i, "/", plot_time_frame, "          ", end="\r")
            timings_str = timings_f.readline()

            if(i>=time_start_frame and time_end_frame):
                timings = np.array(timings_str.split(), dtype=float)
                timings_array[i] = timings[:7]
        labels = ['read_position',\
                  'assemble_rhs',\
                  'precondition',\
                  'mobility_mult',\
                  'k_mult',\
                  'gmres_solve',\
                  'finishing']
        for i in range(len(timings_array[0])):
            ax.plot(time_array, timings_array[:,i], label=labels[i])
            total_timing_array += timings_array[:,i]
        ax.plot(time_array, total_timing_array, label='Total')

        # ax.axvline(299, c='black')
        # plt.annotate('Double precision', (50, 0.055),fontsize=12)
        # plt.annotate('Single precision', (350, 0.035), fontsize=12)
        
        ax.set_ylabel("Computation time/s")
        ax.set_xlabel("Time step")
        ax.set_xlim(time_start_frame, time_end_frame)
        ax.set_ylim(0, 10)
        plt.legend()
        plt.savefig(f'fig/timings_{self.nfil}fil_{self.ar}ar.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()


# Multi sims
    def multi_phase(self):
        # Plotting
        colormap = 'cividis'
        colormap = 'twilight_shifted'

        nrow = len(np.unique(self.pars_list['nfil']))
        ncol = len(np.unique(self.pars_list['ar']))
        # if(ncol == 1 or nrow == 1):
        nrow = int(self.num_sim**.5)
        ncol = nrow + (1 if nrow**2 < self.num_sim else 0)
        print(f'nrow = {nrow} ncol = {ncol}')
        spring_factor = self.pars_list['spring_factor'][0]

        
        fig, axs = plt.subplots(nrow, ncol, figsize=(18, 18), sharex=True, sharey=True)
        # cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # [left, bottom, width, height] for the colorbar

        axs_flat = axs.ravel()
        import scipy.interpolate

        fig.supxlabel(r"Azimuth position")
        fig.supylabel(r"Polar position")
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        # plt.ylabel(r"$\theta$")
        # plt.xlabel(r"$\phi$")
        plt.xlim(-np.pi, np.pi)
        plt.ylim(0, np.pi)
        plt.xticks(np.linspace(-np.pi, np.pi, 5), ['-π', '-π/2', '0', 'π/2', 'π'])
        plt.yticks(np.linspace(0, np.pi, 3), ['0', 'π/2', 'π'])

        for ind, ax in enumerate(axs_flat):
            if (ind < self.num_sim):
                try:
                    self.index = ind
                    self.select_sim()

                    fil_references_sphpolar = np.zeros((self.nfil,3))
                    fil_phases_f = open(self.simName + '_filament_phases.dat', "r")
                    fil_angles_f = open(self.simName + '_filament_shape_rotation_angles.dat', "r")
                    for i in range(self.plot_end_frame):
                        print(" index ", self.index,  " frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                        fil_phases_str = fil_phases_f.readline()
                        # fil_angles_str = fil_angles_f.readline()
                        if(i==self.plot_end_frame-1):
                            fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                            fil_phases = util.box(fil_phases, 2*np.pi)
                            for i in range(self.nfil):
                                fil_references_sphpolar[i] = util.cartesian_to_spherical(self.fil_references[3*i: 3*i+3])
                                
                            if (self.interpolate):
                                n1, n2 = 30, 30
                                azim_grid = np.linspace(-np.pi, np.pi, n1)
                                polar_grid = np.linspace(0, np.pi, n2)
                                xx, yy = np.meshgrid(azim_grid, polar_grid)
                                zz = scipy.interpolate.griddata((fil_references_sphpolar[:,1],fil_references_sphpolar[:,2]), fil_phases, (xx, yy), method='cubic')
                                ax.scatter(xx, yy, c=zz, cmap=colormap)
                            else:
                            # Individual filaments
                                ax.scatter(fil_references_sphpolar[:,1], fil_references_sphpolar[:,2], c=fil_phases, cmap=colormap)

                            # ax.scatter(fil_references_sphpolar[:,1], fil_references_sphpolar[:,2], c=fil_phases, cmap=colormap)
                    # ax.set_ylabel(r"$\theta$")
                    # ax.set_xlabel(r"$\phi$")
                    # ax.set_xlim(-np.pi, np.pi)
                    # ax.set_ylim(0, np.pi)
                    # ax.set_xticks(np.linspace(-np.pi, np.pi, 5), ['-π', '-π/2', '0', 'π/2', 'π'])
                    # ax.set_yticks(np.linspace(0, np.pi, 5), ['0', 'π/4', 'π/2', '3π/4', 'π'])
                    ax.set_title(f"index={self.index} nfil={self.nfil} AR={self.ar}")
                except:
                    print("WARNING: " + self.simName + " not found.")

        plt.tight_layout()
        # plt.savefig(f'fig/ciliate_multi_phase_elst{spring_factor}.png', bbox_inches = 'tight', format='png')
        plt.savefig(f'fig/ciliate_multi_phase_elst{spring_factor}.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def multi_ciliate(self):
        # Plotting
        nrow = int(self.num_sim**.5)
        ncol = nrow + (1 if nrow**2 < self.num_sim else 0)
        fig, axs = plt.subplots(nrow, ncol, figsize=(18, 18), subplot_kw={'projection': '3d'})
        axs_flat = axs.ravel()

        for ind, ax in enumerate(axs_flat):
            ax.set_proj_type('ortho')
            # ax.set_proj_type('persp', 0.05)  # FOV = 157.4 deg
            # ax.view_init(elev=5., azim=45)
            # ax.dist=20
            ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            # ax.axis('off')
            # ax.grid(False)

            # ax.set_xlim(-1000, 1000)
            # ax.set_ylim(-1000, 1000)
            # ax.set_zlim(-1000, 1000)
            if (ind < self.num_sim):
                try:
                    self.index = ind
                    self.select_sim()

                    # Create the sphere data points
                    num_points = 300
                    u = np.linspace(0, 2 * np.pi, num_points)
                    v = np.linspace(0, np.pi, num_points)
                    x = self.radius * np.outer(np.cos(u), np.sin(v))
                    y = self.radius * np.outer(np.sin(u), np.sin(v))
                    z = self.radius * np.outer(np.ones(np.size(u)), np.cos(v))

                    seg_states_f = open(self.simName + '_seg_states.dat', "r")
                    body_states_f = open(self.simName + '_body_states.dat', "r")

                    for i in range(self.plot_end_frame):
                        print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                        body_states_str = body_states_f.readline()
                        seg_states_str = seg_states_f.readline()

                        if(i==self.plot_end_frame-1):
                            body_states = np.array(body_states_str.split()[1:], dtype=float)
                            seg_states = np.array(seg_states_str.split()[1:], dtype=float)

                            for swim in range(self.nswim):
                                body_pos = body_states[7*swim : 7*swim+3]
                                if(self.big_sphere):
                                    ax.plot_surface(x+body_pos[0], y+body_pos[1], z+body_pos[2], color='grey', alpha=0.5)
                                    # ax.plot_surface(x, y, z, color='grey', alpha=0.5)
                                
                                else:
                                    R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                                    for blob in range(int(self.pars['NBLOB'])):
                                        blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                                        ax.scatter(blob_x, blob_y, blob_z)
                                
                                # Robot arm to find segment position (Ignored plane rotation!)
                                for fil in range(self.nfil):
                                    fil_data = np.zeros((self.nseg, 3))
                                    fil_i = int(3*fil*self.nseg)
                                    fil_data[0] = seg_states[fil_i : fil_i+3]

                                    for seg in range(1, self.nseg):
                                        seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                                        fil_data[seg] = seg_pos
                                    ax.plot(fil_data[:,0], fil_data[:,1], fil_data[:,2], c='black', zorder = 100)
                    ax.set_title(f"nfil={self.nfil} AR={self.ar}")
                except:
                    print("WARNING: " + self.simName + " not found.")

        plt.tight_layout()
        plt.savefig(f'fig/ciliate_multi.png', bbox_inches = 'tight', format='png')
        plt.savefig(f'fig/ciliate_multi.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def multi_ciliate_traj(self):
        # Plotting
        nrow = int(self.num_sim**.5)
        ncol = nrow + (1 if nrow**2 < self.num_sim else 0)
        fig, axs = plt.subplots(nrow, ncol, figsize=(18, 18), subplot_kw={'projection': '3d'})
        axs_flat = axs.ravel()

        for ind, ax in enumerate(axs_flat):
            ax.set_proj_type('ortho')
            if (ind < self.num_sim):
                try:
                    self.index = ind
                    self.select_sim()

                    body_states_f = open(self.simName + '_body_states.dat', "r")
                    time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
                    body_pos_array = np.zeros((self.frames, 3))

                    for i in range(self.plot_end_frame):
                        print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                        body_states_str = body_states_f.readline()

                        if(i>=self.plot_start_frame):
                            body_states = np.array(body_states_str.split()[1:], dtype=float)

                            body_pos_array[i-self.plot_start_frame] = body_states[0 : 3]

                    ax.plot(body_pos_array[:,0], body_pos_array[:,1], body_pos_array[:,2])
                    ax.set_title(f"nfil={self.nfil} AR={self.ar}")
                except:
                    print("WARNING: " + self.simName + " not found.")

        plt.tight_layout()
        plt.savefig(f'fig/ciliate_multi_traj.png', bbox_inches = 'tight', format='png')
        plt.savefig(f'fig/ciliate_multi_traj.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def multi_ciliate_speed(self):
         # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for ind in range(self.num_sim):
            try:
                self.index = ind
                self.select_sim()

                body_vels_f = open(self.simName + '_body_vels.dat', "r")

                time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
        
                body_vel_array = np.zeros((self.frames, 6))
                body_speed_array = np.zeros(self.frames)

                for i in range(self.plot_end_frame):
                    print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                    body_vels_str = body_vels_f.readline()

                    if(i>=self.plot_start_frame):
                        body_vels = np.array(body_vels_str.split(), dtype=float)

                        body_vel_array[i-self.plot_start_frame] = body_vels
                        body_speed_array[i-self.plot_start_frame] = np.sqrt(np.sum(body_vels[0:3]*body_vels[0:3], 0))

                ax.plot(time_array, body_speed_array, label=f"nfil={self.nfil} AR={self.ar}")
            except:
                print("WARNING: " + self.simName + " not found.")

        plt.tight_layout()
        plt.savefig(f'fig/multi_speed.png', bbox_inches = 'tight', format='png')
        plt.savefig(f'fig/multi_speed.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def multi_timing(self):
        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for ind in range(self.num_sim):
            try:
                self.index = ind
                self.select_sim()

                with open(self.simName + '_time.dat', 'r') as file:
                    plot_time_frame = min(300, len(file.readlines()))
                timings_f = open(self.simName + '_time.dat', "r")

                time_array = np.arange(0, plot_time_frame)
                timings_array = np.zeros((len(time_array), 7))
                total_timing_array = np.zeros(len(time_array))

                for i in range(plot_time_frame):
                    print(" frame ", i, "/", plot_time_frame, "          ", end="\r")
                    timings_str = timings_f.readline()

                    timings = np.array(timings_str.split(), dtype=float)

                    timings_array[i] = timings[:7]
                for i in range(len(timings_array[0])):
                    # ax.plot(time_array, timings_array[:,i], label=labels[i])
                    total_timing_array += timings_array[:,i]
                ax.plot(time_array, total_timing_array, label=f"nfil={self.nfil} AR={self.ar}")
            except:
                print("WARNING: " + self.simName + " not found.")

        plt.tight_layout()
        plt.savefig(f'fig/multi_timing.png', bbox_inches = 'tight', format='png')
        plt.savefig(f'fig/multi_timing.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def multi_ciliate_svd(self):
         # Plotting
        nrow = len(np.unique(self.pars_list['nfil']))
        ncol = len(np.unique(self.pars_list['ar']))
        spring_factor = self.pars_list['spring_factor'][0]

        # nrow = int(self.num_sim**.5)
        # ncol = nrow + (1 if nrow**2 < self.num_sim else 0)
        fig, axs = plt.subplots(nrow, ncol, figsize=(18, 18), sharex=True, sharey=True)
    
        axs_flat = axs.ravel()

        for ind, ax in enumerate(axs_flat):
            if (ind < self.num_sim):
                try:
                    self.index = ind
                    self.select_sim()

                    fil_phases_f = open(self.simName + '_filament_phases.dat', "r")
                    fil_angles_f = open(self.simName + '_filament_shape_rotation_angles.dat', "r")


                    nfil = self.nfil
                    n_snapshots = min(900, self.plot_end_frame)
                    start = self.plot_end_frame - n_snapshots
                    X = np.zeros((nfil, n_snapshots))
                    X_angle = np.zeros((nfil, n_snapshots))

                    fil_references_sphpolar = np.zeros((nfil,3))
                    for fil in range(nfil):
                        fil_references_sphpolar[fil] = util.cartesian_to_spherical(self.fil_references[3*fil: 3*fil+3])
                    azim_array = fil_references_sphpolar[:,1]
                    polar_array = fil_references_sphpolar[:,2]
                    sorted_indices = np.argsort(azim_array)
                    azim_array_sorted = azim_array[sorted_indices]
                    polar_array_sorted = polar_array[sorted_indices]

                    for i in range(self.plot_end_frame):
                        print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                        fil_phases_str = fil_phases_f.readline()
                        fil_angles_str = fil_angles_f.readline()
                        
                        if(i>=start):
                            fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                            fil_phases_sorted = fil_phases[sorted_indices]
                            fil_phases_sorted = util.box(fil_phases_sorted, 2*np.pi)

                            fil_phases_sorted = np.sin(fil_phases_sorted)

                            fil_angles = np.array(fil_angles_str.split()[1:], dtype=float)
                            fil_angles_sorted = fil_angles[sorted_indices]

                            X[:,i-start] = fil_phases_sorted[:nfil]
                            X_angle[:,i-start] = fil_angles_sorted[:nfil]
                    
                    np.savetxt(f'phase_data/by_index/spring_constant{spring_factor}/X_phase_index{self.index}.txt', X, delimiter=', ')
                    np.savetxt(f'phase_data/by_index/spring_constant{spring_factor}/X_rotation_angle_index{self.index}.txt', X_angle, delimiter=', ')
                    np.savetxt(f'phase_data/by_index/spring_constant{spring_factor}/azim_pos_index{self.index}.txt', azim_array_sorted, delimiter=', ')
                    np.savetxt(f'phase_data/by_index/spring_constant{spring_factor}/polar_pos_index{self.index}.txt', polar_array_sorted, delimiter=', ')

                    np.savetxt(f'phase_data/by_pars/spring_constant{spring_factor}/X_phase_nfil{nfil}_rol{self.ar}_spring{self.spring_factor}.txt', X, delimiter=', ')
                    np.savetxt(f'phase_data/by_pars/spring_constant{spring_factor}/X_rotation_angle_nfil{nfil}_rol{self.ar}_spring{self.spring_factor}.txt', X_angle, delimiter=', ')
                    np.savetxt(f'phase_data/by_pars/spring_constant{spring_factor}/azim_pos_nfil{nfil}_rol{self.ar}_spring{self.spring_factor}.txt', azim_array_sorted, delimiter=', ')
                    np.savetxt(f'phase_data/by_pars/spring_constant{spring_factor}/polar_pos_nfil{nfil}_rol{self.ar}_spring{self.spring_factor}.txt', polar_array_sorted, delimiter=', ')


                    U, sigma, V = np.linalg.svd(X, full_matrices=False)
                    Sigma = np.diag(sigma)

                    pc = U @ Sigma
                    pa = V

                    # for fil in range(num_fil):
                    #     abs_pc = np.abs(pc[fil][:nfil])
                    #     ax.plot(np.cumsum(abs_pc)/np.sum(abs_pc), label=f'fil {fil}')
                    # ax.set_xlabel('Mode')
                    # ax.set_ylabel('Accumulated |weight| fraction')
                    # ax.legend()

                except:
                    print("WARNING: " + self.simName + " not found.")

        # plt.tight_layout()
        # plt.savefig(f'fig/multi_svd.png', bbox_inches = 'tight', format='png')
        # plt.savefig(f'fig/multi_svd.pdf', bbox_inches = 'tight', format='pdf')
        # plt.show()
    
   

# Summary plot
    def summary_ciliate_speed(self):
        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1,1,1)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(1,1,1)

        nfil_list = np.array(self.pars_list['nfil'])
        ar_list = np.array(self.pars_list['ar'])
        sphere_r_list = np.zeros(self.num_sim)
        speed_list = np.zeros(self.num_sim)
        angular_speed_list = np.zeros(self.num_sim)

        for ind in range(self.num_sim):
            try:
                self.index = ind
                self.select_sim()

                body_vels_f = open(self.simName + '_body_vels.dat', "r")

                body_speed_array = np.zeros(self.frames)
                body_angular_speed_array = np.zeros(self.frames)

                for i in range(self.plot_end_frame):
                    print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                    body_vels_str = body_vels_f.readline()

                    if(i>=self.plot_start_frame):
                        body_vels = np.array(body_vels_str.split(), dtype=float)

                        body_speed_array[i-self.plot_start_frame] = np.sqrt(np.sum(body_vels[0:3]*body_vels[0:3], 0))
                        body_angular_speed_array[i-self.plot_start_frame] = np.sqrt(np.sum(body_vels[3:6]*body_vels[3:6], 0))

                speed_list[ind] = np.mean(body_speed_array)
                angular_speed_list[ind] = np.mean(body_angular_speed_array)
                sphere_r_list[ind] = self.radius
            except:
                print("WARNING: " + self.simName + " not found.")
        
        colormap = 'Greys'
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        vmin, vmax = np.min(speed_list), np.max(speed_list)
        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig2.colorbar(sm)
        cbar.ax.set_yticks(np.linspace(vmin, vmax, 8))
        cbar.set_label(r"Speed")

        ax2.scatter(self.pars_list['nfil'], self.pars_list['ar'], s=100, c=speed_list, edgecolors='red', linewidths=0.1, cmap=colormap)
        ax2.set_xlabel("Nfil")
        ax2.set_ylabel("R/L")
        fig2.savefig(f'fig/multi_ciliate_speed_summary_heatmap.png', bbox_inches = 'tight', format='png')
        fig2.savefig(f'fig/multi_ciliate_speed_summary_heatmap.pdf', bbox_inches = 'tight', format='pdf')

        vmin, vmax = np.min(angular_speed_list), np.max(angular_speed_list)
        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig3.colorbar(sm)
        cbar.ax.set_yticks(np.linspace(vmin, vmax, 8))
        cbar.set_label(r"Angular speed")

        ax3.scatter(self.pars_list['nfil'], self.pars_list['ar'], s=100, c=angular_speed_list, edgecolors='red', linewidths=0.1, cmap=colormap)
        ax3.set_xlabel("Nfil")
        ax3.set_ylabel("R/L")
        fig3.savefig(f'fig/multi_ciliate_angular_speed_summary_heatmap.png', bbox_inches = 'tight', format='png')
        fig3.savefig(f'fig/multi_ciliate_angular_speed_summary_heatmap.pdf', bbox_inches = 'tight', format='pdf')
        

        ax.scatter(nfil_list, speed_list, label=f"nfil={self.nfil} AR={self.ar}")
        ax.set_ylabel("Velocity")
        ax.set_xlabel("Number of filaments")
        fig.tight_layout()
        fig.savefig(f'fig/multi_ciliate_speed_summary.png', bbox_inches = 'tight', format='png')
        fig.savefig(f'fig/multi_ciliate_speed_summary.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def summary_timing(self):
        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        num_particle_list = np.array(self.pars_list['nfil'])*20 + np.array(self.pars_list['nblob'])
        total_timing_list = np.zeros((self.num_sim, 7))

        for ind in range(self.num_sim):
            try:
                self.index = ind
                self.select_sim()

                with open(self.simName + '_time.dat', 'r') as file:
                    plot_time_frame = min(300, len(file.readlines()))
                timings_f = open(self.simName + '_time.dat', "r")

                for i in range(plot_time_frame):
                    print(" frame ", i, "/", plot_time_frame, "          ", end="\r")
                    timings_str = timings_f.readline()

                    timings = np.array(timings_str.split()[:-1], dtype=float)

                    total_timing_list[ind] += timings
            except:
                print("WARNING: " + self.simName + " not found.")


        ax.scatter(num_particle_list, total_timing_list[:,3]/plot_time_frame, label="HI solver")
        ax.scatter(num_particle_list, np.sum(total_timing_list, axis=1)/plot_time_frame, label="Total")
        # ax.set_ylim(0)
        # ax.set_yscale('log')
        ax.set_ylabel("Compute time per time step/s")
        ax.set_xlabel("Number of particles")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'fig/multi_timing_summary.png', bbox_inches = 'tight', format='png')
        plt.savefig(f'fig/multi_timing_summary.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def summary_ciliate_dissipation(self):
        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1,1,1)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(1,1,1)
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(1,1,1)
        fig5 = plt.figure()
        ax5 = fig5.add_subplot(1,1,1)

        nfil_list = np.array(self.pars_list['nfil'])
        ar_list = np.array(self.pars_list['ar'])
        sphere_r_list = np.zeros(self.num_sim)
        speed_list = np.zeros(self.num_sim)
        angular_speed_list = np.zeros(self.num_sim)
        dissipation_list = np.zeros(self.num_sim)
        efficiency_list = np.zeros(self.num_sim)

        for ind in range(self.num_sim):
            try:
                self.index = ind
                self.select_sim()

                seg_forces_f = open(self.simName + '_seg_forces.dat', "r")
                seg_vels_f = open(self.simName + '_seg_vels.dat', "r")
                blob_forces_f = open(self.simName + '_blob_forces.dat', "r")
                blob_references_f = open(self.simName + '_blob_references.dat', "r")
                body_vels_f = open(self.simName + '_body_vels.dat', "r")
                
                body_speed_array = np.zeros(self.frames)
                body_angular_speed_array = np.zeros(self.frames)
                dissipation_array = np.zeros(self.frames)

                blob_references_str = blob_references_f.readline()
                blob_references = np.array(blob_references_str.split(), dtype=float)
                blob_references = np.reshape(blob_references, (int(self.pars['NBLOB']), 3))
            
                for i in range(self.plot_end_frame):
                    print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                    seg_forces_str = seg_forces_f.readline()
                    seg_vels_str = seg_vels_f.readline()
                    blob_forces_str = blob_forces_f.readline()
                    body_vels_str = body_vels_f.readline()

                    if(i>=self.plot_start_frame):
                        seg_forces = np.array(seg_forces_str.split()[1:], dtype=float)
                        seg_vels = np.array(seg_vels_str.split()[1:], dtype=float)
                        blob_forces= np.array(blob_forces_str.split()[1:], dtype=float)
                        body_vels= np.array(body_vels_str.split(), dtype=float)

                        seg_forces = np.reshape(seg_forces, (int(self.pars['NSEG']*self.pars['NFIL']), 6))
                        seg_vels = np.reshape(seg_vels, (int(self.pars['NSEG']*self.pars['NFIL']), 6))
                        blob_forces = np.reshape(blob_forces, (int(self.pars['NBLOB']), 3))
                        body_vels_tile = np.tile(body_vels, (int(self.pars['NBLOB']), 1))
                        blob_vels = body_vels_tile[:, 0:3] + np.cross(body_vels_tile[:, 3:6], blob_references)

                        body_speed_array[i-self.plot_start_frame] = np.sqrt(np.sum(body_vels[0:3]*body_vels[0:3], 0))
                        body_angular_speed_array[i-self.plot_start_frame] = np.sqrt(np.sum(body_vels[3:6]*body_vels[3:6], 0))
                        dissipation_array[i-self.plot_start_frame] = np.sum(blob_forces * blob_vels) + np.sum(seg_forces * seg_vels)

                speed_list[ind] = np.mean(body_speed_array)
                angular_speed_list[ind] = np.mean(body_angular_speed_array)
                dissipation_list[ind] = np.mean(dissipation_array)
                sphere_r_list[ind] = self.radius
            except:
                print("WARNING: " + self.simName + " not found.")

        efficiency_list = 6*np.pi*sphere_r_list*speed_list**2/dissipation_list
        print(repr(speed_list))
        print(repr(dissipation_list))
        print(repr(efficiency_list))

        colormap = 'Greys'
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable

        # Speed
        vmin, vmax = np.min(speed_list), np.max(speed_list)
        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig2.colorbar(sm)
        cbar.ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
        cbar.ax.set_yticks(np.linspace(vmin, vmax, 8))
        cbar.set_label(r"Speed")

        ax2.scatter(self.pars_list['nfil'], self.pars_list['ar'], s=100, c=speed_list, edgecolors='red', linewidths=0.1, cmap=colormap)
        ax2.set_xlabel("Nfil")
        ax2.set_ylabel("R/L")
        # fig2.savefig(f'fig/multi_ciliate_speed_summary_heatmap.png', bbox_inches = 'tight', format='png')
        fig2.savefig(f'fig/multi_ciliate_speed_summary_heatmap.pdf', bbox_inches = 'tight', format='pdf')

        # Angular speed
        vmin, vmax = np.min(angular_speed_list), np.max(angular_speed_list)
        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig3.colorbar(sm)
        cbar.ax.set_yticks(np.linspace(vmin, vmax, 8))
        cbar.set_label(r"Angular speed")

        ax3.scatter(self.pars_list['nfil'], self.pars_list['ar'], s=100, c=angular_speed_list, edgecolors='red', linewidths=0.1, cmap=colormap)
        ax3.set_xlabel("Nfil")
        ax3.set_ylabel("R/L")
        # fig3.savefig(f'fig/multi_ciliate_angular_speed_summary_heatmap.png', bbox_inches = 'tight', format='png')
        fig3.savefig(f'fig/multi_ciliate_angular_speed_summary_heatmap.pdf', bbox_inches = 'tight', format='pdf')

        # Dissipation
        vmin, vmax = np.min(dissipation_list), np.max(dissipation_list)
        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig4.colorbar(sm)
        cbar.ax.set_yticks(np.linspace(vmin, vmax, 8))
        cbar.set_label(r"Dissipation")

        ax4.scatter(self.pars_list['nfil'], self.pars_list['ar'], s=100, c=dissipation_list, edgecolors='red', linewidths=0.1, cmap=colormap)
        ax4.set_xlabel("Nfil")
        ax4.set_ylabel("R/L")
        # fig4.savefig(f'fig/multi_ciliate_dissipation_summary_heatmap.png', bbox_inches = 'tight', format='png')
        fig4.savefig(f'fig/multi_ciliate_dissipation_summary_heatmap.pdf', bbox_inches = 'tight', format='pdf')

        # Efficiency
        vmin, vmax = np.min(efficiency_list), np.max(efficiency_list)
        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig5.colorbar(sm)
        cbar.ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
        cbar.ax.set_yticks(np.linspace(vmin, vmax, 8))
        cbar.set_label(r"Efficiency")

        ax5.scatter(self.pars_list['nfil'], self.pars_list['ar'], s=100, c=efficiency_list, edgecolors='red', linewidths=0.1, cmap=colormap)
        ax5.set_xlabel("Nfil")
        ax5.set_ylabel("R/L")
        # fig5.savefig(f'fig/multi_ciliate_efficiency_summary_heatmap.png', bbox_inches = 'tight', format='png')
        fig5.savefig(f'fig/multi_ciliate_efficiency_summary_heatmap.pdf', bbox_inches = 'tight', format='pdf')
        


        ax.scatter(nfil_list, dissipation_list, label=f"nfil={self.nfil} AR={self.ar}")
        ax.set_ylabel("Dissipation")
        ax.set_xlabel("Number of filaments")
        fig.tight_layout()
        # fig.savefig(f'fig/ciliate_dissipation_summary.png', bbox_inches = 'tight', format='png')
        fig.savefig(f'fig/ciliate_dissipation_summary.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def summary_check_overlap(self):
        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for ind in range(self.num_sim):
            # try:
            self.index = ind
            self.select_sim()

            seg_states_f = open(self.simName + '_seg_states.dat', "r")
            body_states_f = open(self.simName + '_body_states.dat', "r")

            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.plot_end_frame, "          ", end="\r")
                segs_list = np.zeros((int(self.nfil*self.nseg), 3))
                blobs_list = np.zeros((int(self.nblob), 3))

                body_states_str = body_states_f.readline()
                seg_states_str = seg_states_f.readline()
                
                if(i==self.plot_end_frame-1):
                    body_states = np.array(body_states_str.split()[1:], dtype=float)
                    seg_states = np.array(seg_states_str.split()[1:], dtype=float)

                    for swim in range(int(self.pars['NSWIM'])):
                        body_pos = body_states[7*swim : 7*swim+3]
                        R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                        R = np.linalg.inv(R)
                        R = np.eye(3)
                        for blob in range(int(self.pars['NBLOB'])):
                            blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                            blobs_list[blob] = blob_x, blob_y, blob_z
                        for fil in range(int(self.pars['NFIL'])):
                            fil_i = int(3*fil*self.pars['NSEG'])
                            old_seg_pos = seg_states[fil_i : fil_i+3]
                            segs_list[fil*self.nseg] = old_seg_pos
                            for seg in range(1, int(self.pars['NSEG'])):
                                seg_pos = seg_states[fil_i+3*(seg) : fil_i+3*(seg+1)]
                                segs_list[fil*self.nseg + seg] = seg_pos
                    threshold = 0.95
                    colliding_indices, colliding_particles = util.label_colliding_particles_with_3d_cell_list(segs_list, 5, threshold*float(self.pars['RSEG']))
                    if not colliding_indices:
                        print(f'index={ind} - No overlapping at threshold {threshold}\n')
                    else:
                        print(f'index={ind} - Overlapping at threshold {threshold}\n')
            # except:
            #     print("WARNING: " + self.simName + " not found.")
        
        # ax.scatter(nfil_list, speed_list, label=f"nfil={self.nfil} AR={self.ar}")
        # ax.set_ylabel("Velocity")
        # ax.set_xlabel("Number of filaments")
        # fig.tight_layout()
        # fig.savefig(f'fig/multi_ciliate_speed_summary.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()
# 
# density for different elsticity
# efficiency/disspation
# one big sim start from random
# phi dot data - constant omega
# 

# Special plot
    def ishikawa(self):
        top_dir = "data/expr_sims/ishikawa_expr/"
        vel_dirs = dirs = [ "k0.0/",
                            "k0.5/",
                            "k1.0/",
                            "k1.5/",
                            "k2.0/"]
        
        # dissipation_dirs = ["k0.0/",
        #                     "k1.0/"]
        
        N_dirs = ["k0.0N162/",
                  "k0.0N636/",
                  "k0.0N2520/",]
        # N_dirs = ["k0.0N162_RoL20/",
        #           "k0.0N636_RoL20/",
        #           "k0.0N2520_RoL20/",]

        dirs = ["k0.0/",
                "k0.5/",
                "k1.0/",
                "k1.5/",
                "k2.0/",
                "k0.0N162/",
                "k0.0N636/",
                "k0.0N2520/"]
        
        ls = ['solid', 'dashed', 'dotted']
        markers = ["^", "s", "d"]
        labels = [r"$k=0$",r"$k=0.5$",r"$k=1$",r"$k=1.5$",r"$k=2$",]
        colors = ["black","red","green","blue","purple"]
        dissipation_labels = [r"$N=162$",r"$N=636$",r"$N=2520$"]
        dissipation_colors = ["green","black","red"]
        vel_marker = 0
        dissipation_marker = 0

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_xlim(0, 1)
        ax.set_ylim(-1.4, 4.5)
        ax.set_ylabel(r"$V_zT/L$")
        ax.set_xlabel(r"$t/T$")

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1,1,1)
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 80000)
        ax2.set_ylabel(r"$PT^2/\mu L^3$")
        ax2.set_xlabel(r"$t/T$")

        for ni, directory in enumerate(N_dirs):
            try:
                self.dir = top_dir + directory
                self.read_rules()
                self.select_sim()
                L = self.pars['FIL_LENGTH']
                print(f'dt={self.dt}')
                print(f'L={L}')
                start_time = 0
                end_time = min(32, sum(1 for line in open(self.simName + '_body_states.dat')))

                body_states_f = open(self.simName + '_body_states.dat', "r")
                body_vels_f = open(self.simName + '_body_vels.dat', "r")
                seg_forces_f = open(self.simName + '_seg_forces.dat', "r")
                seg_vels_f = open(self.simName + '_seg_vels.dat', "r")
                blob_forces_f = open(self.simName + '_blob_forces.dat', "r")
                blob_references_f = open(self.simName + '_blob_references.dat', "r")

                time_array = np.arange(start_time, end_time)
                body_pos_array = np.zeros((len(time_array), 3))
                body_vel_array = np.zeros((len(time_array), 6))
                body_speed_array = np.zeros(len(time_array))
                body_angular_speed_array = np.zeros(len(time_array))
                dissipation_array = np.zeros(len(time_array))

                blob_references_str = blob_references_f.readline()
                blob_references = np.array(blob_references_str.split(), dtype=float)
                blob_references = np.reshape(blob_references, (int(self.pars['NBLOB']), 3))

                for i in range(end_time):
                    print(" frame ", i, "/", end_time, "          ", end="\r")
                    seg_forces_str = seg_forces_f.readline()
                    seg_vels_str = seg_vels_f.readline()
                    blob_forces_str = blob_forces_f.readline()
                    body_vels_str = body_vels_f.readline()

                    if(i>=start_time):
                        seg_forces = np.array(seg_forces_str.split()[1:], dtype=float)
                        seg_vels = np.array(seg_vels_str.split()[1:], dtype=float)
                        blob_forces= np.array(blob_forces_str.split()[1:], dtype=float)
                        body_vels = np.array(body_vels_str.split(), dtype=float)
                        

                        seg_forces = np.reshape(seg_forces, (int(self.pars['NSEG']*self.pars['NFIL']), 6))
                        seg_vels = np.reshape(seg_vels, (int(self.pars['NSEG']*self.pars['NFIL']), 6))
                        blob_forces = np.reshape(blob_forces, (int(self.pars['NBLOB']), 3))
                        body_vels_tile = np.tile(body_vels, (int(self.pars['NBLOB']), 1))
                        blob_vels = body_vels_tile[:, 0:3] + np.cross(body_vels_tile[:, 3:6], blob_references)

                        body_vel_array[i-start_time] = body_vels
                        body_speed_array[i-start_time] = np.sqrt(np.sum(body_vels[0:3]*body_vels[0:3], 0))
                        body_angular_speed_array[i-start_time] = np.sqrt(np.sum(body_vels[3:6]*body_vels[3:6], 0))
                        dissipation_array[i-start_time] = np.sum(blob_forces * blob_vels) + np.sum(seg_forces[:, 0:3] * seg_vels[:, 0:3])

                # if directory in vel_dirs:
                #     ax.plot(time_array/30., body_vel_array[:,2]/L, label=labels[vel_marker], c=colors[vel_marker])
                #     vel_marker += 1
                if directory in N_dirs:
                    ax.plot(time_array*self.dt, body_vel_array[:,2]/L, label=dissipation_labels[dissipation_marker], c=dissipation_colors[dissipation_marker])
                    ax2.plot(time_array*self.dt, dissipation_array/L**3, label=dissipation_labels[dissipation_marker], c=dissipation_colors[dissipation_marker])
                    dissipation_marker +=1
            except:
                print("WARNING: " + self.dir + "rules.ini not found.")
        
        # Plot the comparison data
        # JFM 2019 velocity
        directory = 'pyfile/analysis/ishikawa_data/'
        files = ['k0.0.csv', 'k0.5.csv', 'k1.0.csv', 'k1.5.csv', 'k2.0.csv']
        files = ['vel_k0.0N162.csv', 'vel_k0.0N636.csv', 'vel_k0.0N2520.csv']
        for i, filename in enumerate(files):
            try:
                file = open(directory + filename, mode='r')
                df = pd.read_csv(directory + filename, header=None)
                data = df.to_numpy()
                x, y = data[:,0], data[:,1]

                ax.plot(x, y, ls='dotted', c=dissipation_colors[i], alpha=0.5)
            except:
                print("WARNING: " + directory + filename + " not found.")
        # PNAS 2020 dissipation
        files = ['dissipation_k0.0N162.csv', 'dissipation_k0.0N636.csv', 'dissipation_k0.0N2520.csv']
        for i, filename in enumerate(files):
            try:
                file = open(directory + filename, mode='r')
                df = pd.read_csv(directory + filename, header=None)
                data = df.to_numpy()
                x, y = data[:,0], data[:,1]

                ax2.plot(x, y, ls='dotted', c=dissipation_colors[i], alpha=0.5)
            except:
                print("WARNING: " + directory + filename + " not found.")


        # Make legends
        legend_label = r'$Ito\ etc.\ (2019)$'
        legend_label = r'$Omori\ etc.\ (2020)$'
        legend1 = ax.legend()
        line1a, = ax.plot([-1, -1.1], [-1, -1.1], ls='-', c='black', label=r'$data$' )
        line1b, = ax.plot([-1, -1.1], [-1, -1.1], ls='dotted', c='black', label=legend_label)
        ax.legend(handles = [line1a, line1b], loc='upper left')
        ax.add_artist(legend1)
        
        legend2 = ax2.legend()
        line2a, = ax2.plot([-1, -1.1], [-1, -1.1], ls='-', c='black', label=r'$data$' )
        line2b, = ax2.plot([-1, -1.1], [-1, -1.1], ls='dotted', c='black', label=legend_label)
        ax2.legend(handles = [line2a, line2b], loc='upper left')
        ax2.add_artist(legend2)

        fig.savefig(f'fig/PNAS_comparison_vel.pdf', bbox_inches = 'tight', format='pdf')
        fig2.savefig(f'fig/PNAS_comparison_dissipation.pdf', bbox_inches = 'tight', format='pdf')
        # fig.savefig(f'fig/ishikawa_vel{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        # fig2.savefig(f'fig/ishikawa_dissipation{self.nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

#