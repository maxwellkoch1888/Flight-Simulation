import json
import numpy as np
import matplotlib.pyplot as plt
import time 

from quaternion import quat_mult, quat_norm, quat_to_euler, eul_to_quat

# OPEN THE JSON FILE
with open('json/graphics.json', 'r') as file:
    graphics = json.load(file)

# UNPACK THE VALUES IN JSON FILE  
dvp   = graphics['camera']['view_plane']['distance[ft]']
theta = graphics['camera']['view_plane']['angle[deg]']
RAvp  = graphics['camera']['view_plane']['aspect_ratio']

camera_location    = graphics['camera']['location[ft]']
camera_orientation = graphics['camera']['orientation[deg]']

altitude = graphics["ground"]["altitude[ft]"]
grid_number = graphics["ground"]["grid_number"]
grid_scale = graphics["ground"]["grid_scale[ft]"]
color = graphics["ground"]["color"]


camera_location = np.array([[camera_location[0]], [camera_location[1]], [camera_location[2]]])
theta              = np.radians(theta)
camera_orientation = np.radians(camera_orientation)
camera_orientation_quat = eul_to_quat(camera_orientation)


class Camera:
    def __init__(self, orientation:np.ndarray, location:np.ndarray, vp_distance:float, vp_angle:float, vp_aspect_ratio:float):
        # FIND THE CAMERA LOCATION AND ORIENTATION
        eul = orientation
        self.location = location 
        self.quat = eul_to_quat(eul)

        # FIND THE VIEW PLANE PROPERTIES
        self.vp_distance = vp_distance
        self.vp_angle = vp_angle
        self.vp_aspect_ratio = vp_aspect_ratio

        self.wvp = 2.0 * self.vp_distance * np.tan(self.vp_angle)
        self.hvp = self.wvp / self.vp_aspect_ratio

        # INITIALIZE VIEW PLANE BODY FIXED COORDINATES
        self.vp_xv = np.array([self.vp_distance, self.vp_distance, self.vp_distance, self.vp_distance], dtype=np.float64)
        self.vp_yv = np.array([-self.wvp, -self.wvp, self.wvp, self.wvp], dtype=np.float64) * 0.5
        self.vp_zv = np.array([-self.hvp, self.hvp, self.hvp, -self.hvp], dtype=np.float64) * 0.5 
        camera_fixed = np.array([self.vp_xv, self.vp_yv, self.vp_zv])

        # INITIALIZE VIEW PLANE EARTH FIXED COORDINATES 
        self.vp_xf = np.array([0.0, 0.0, 0.0, 0.0])
        self.vp_yf = np.array([0.0, 0.0, 0.0, 0.0])
        self.vp_zf = np.array([0.0, 0.0, 0.0, 0.0])
        
        self.lines2D = np.zeros((2,2))

        self.dx = self.wvp
        self.dy = self.hvp

    def set_state(self, location:np.ndarray, quat:np.ndarray):
        # SET THE LOCATION AND ORIENTATION 
        self.location = location 
        self.quat = quat 

        # CALCULATE VIEW IN EARTH FIXED PLANE
        e0 = self.quat[0]
        ex = self.quat[1]
        ey = self.quat[2]
        ez = self.quat[3]

        camera_fixed = np.array([self.vp_xv, self.vp_yv, self.vp_zv])

        rotation_matrix = np.array([[ex**2 + e0**2 - ey**2 - ez**2, 2*(ex*ey + ez*e0), 2*(ex*ez - ey*e0)],
                                    [2*(ex*ey - ez*e0), ey**2 + e0**2 - ex**2 - ez**2, 2*(ey*ez + ex*e0)],
                                    [2*(ex*ez + ey*e0), 2*(ey*ez - ex*e0), ez**2 + e0**2 - ex**2 - ey**2]])

        earth_fixed = np.matmul(np.linalg.inv(rotation_matrix), camera_fixed) + camera_location            

        self.vp_xf = earth_fixed[0,:]
        self.vp_yf = earth_fixed[1,:]
        self.vp_zf = earth_fixed[2,:]

        # BUILD CENTER POINT OF PLANE
        self.P0 = np.array([[np.mean(earth_fixed[0])], [np.mean(earth_fixed[1])], [np.mean(earth_fixed[2])]]).flatten()

        # CALCULATE THE NORMAL VECTOR 
        P1 = np.array([earth_fixed[0,0], earth_fixed[1,0], earth_fixed[2,0]])
        P2 = np.array([earth_fixed[0,1], earth_fixed[1,1], earth_fixed[2,1]])
        self.nvp = np.cross(P1 - self.P0, P2 - self.P0)
            
class Lines:
    def __init__(self, grid_number:int, grid_scale:float, altitude:float, ax):
        # BUILD POINTS TO BE PLOTTED 
        nxlines = grid_number * 2 + 1 
        self.nlines = 2*nxlines
        self.points = np.zeros((4*nxlines,3))
        self.lines = np.zeros((2*nxlines,2), dtype=int)
        for i in range(nxlines):
            self.points[2*i,   :] = [-grid_number*grid_scale, (i-grid_number)*grid_scale, -altitude]
            self.points[2*i+1, :] = [ grid_number*grid_scale, (i-grid_number)*grid_scale, -altitude]
            self.lines[i,0] = 2*i
            self.lines[i,1] = 2*i+1

        for i in range(nxlines):
            self.points[2*i+2*nxlines  , :] = [(i-grid_number)*grid_scale, -grid_number*grid_scale, -altitude]
            self.points[2*i+2*nxlines+1, :] = [(i-grid_number)*grid_scale,  grid_number*grid_scale, -altitude]
            self.lines[i+nxlines,0] = 2*i+2*nxlines 
            self.lines[i+nxlines,1] = 2*i+2*nxlines+1

        self.npoints = len(self.points)
        self.nlines = len(self.lines)

        # SET UP AXES FOR DRAWING
        self.ax, = ax.plot([], [], marker='', ls='-', color=color)

        # ALLOCATE SPACE FOR 3D POINTS IN EARTH FIXED COORDINATES
        self.points3D = np.zeros((self.npoints,3))
        # self.set_state([0.0, 0.0, 0.0,], [1.0, 0.0, 0.0, 0.0])

        # ALLOCATE SPACE FOR 2D POINTS ON VIEWPLANE
        self.points2D = np.zeros((self.npoints,2))
        self.lca = np.zeros((self.npoints,3))
        self.t = np.zeros(self.npoints)
        self.pf = np.zeros((self.npoints,3))

        self.lines2D = np.full((self.nlines*3, 2), np.nan, dtype=float)

        self.points3D[:] = self.points

    def plot(self, camera):
        pc = camera.location
        lc0 = (camera.P0 - pc.flatten()).flatten()
        temp = np.dot(lc0, camera.nvp)

        self.lca[:] = self.points3D[:] - pc.flatten()

        self.t[:] = temp/np.dot(self.lca,camera.nvp)

        self.pf[:] = self.lca[:]*self.t[:, None]

        e0, ex, ey, ez = camera.quat
        a00 = ex*ex + e0*e0 - ey*ey - ez*ez 
        a01 = 2*(ex*ey - ez*e0)
        a02 = 2*(ex*ez - ey*e0)
        a10 = 2*(ex*ey - ez*e0)
        a11 = ey*ey + e0*e0 - ex*ex - ez*ez 
        a12 = 2*(ey*ez - ex*e0)
        a20 = 2*(ex*ez - ey*e0)
        a21 = 2*(ey*ez - ex*e0)
        a22 = ez*ez + e0*e0 - ey*ey - ex*ex

        self.points2D[:,0] =  a01*self.pf[:,0] + a11*self.pf[:,1] + a21*self.pf[:,2]
        self.points2D[:,1] = -a02*self.pf[:,0] - a12*self.pf[:,1] - a22*self.pf[:,2]
        
        for i in range(self.nlines):
            i0 = self.lines[i,0]
            i1 = self.lines[i,1]

            if(self.t[i0]>0 and self.t[i1]>0):
                self.lines2D[3*i,  :] = self.points2D[i0,:]
                self.lines2D[3*i+1,:] = self.points2D[i1,:]
            # elif(self.t[i0]<0 and self.t[i1]<0):

        self.ax.set_data(self.lines2D[:,0], self.lines2D[:,1])

# BUILD THE CAMERA
camera = Camera(
            orientation= camera_orientation,
            location=camera_location, 
            vp_distance=dvp,
            vp_angle=theta,
            vp_aspect_ratio=RAvp
            )

quat_init = eul_to_quat(camera_orientation)

camera.set_state(
    location=camera_location, 
    quat=quat_init
    )

# PLOT THE FIGURE
fig = plt.figure(figsize=(camera.vp_aspect_ratio*5.0,5.0))
ax = fig.add_subplot(111)
plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
plt.axis('off')
ax.axes.set_xlim(camera.vp_yv[0], camera.vp_yv[2])
ax.axes.set_ylim(-camera.vp_zv[1], -camera.vp_zv[0])
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.set_xticks([])
ax.set_yticks([])
ax.axes.set_aspect('equal')
fig.canvas.draw()
# ground.plot()
# plt.show()
plt.show(block=False)

# BUILD THE LINES
ground = Lines(
            grid_number=grid_number,
            grid_scale=grid_scale,
            altitude=altitude,
            ax=ax
            )


# PRINT RESULTS
print('------------------------------------ Results ------------------------------------')
print('11.2.1:')
print(f'w_vp:    {camera.wvp:.11f} ft')
print(f'h_vp:    {camera.hvp:.11f} ft')
print('')
np.set_printoptions(precision=2, suppress=True)
print(f'x_c_vp:  {camera.vp_xv} ft')
np.set_printoptions(precision=11, suppress=True)
print(f'y_c_vp:  {camera.vp_yv} ft')
np.set_printoptions(precision=12, suppress=True)
print(f'z_c_vp:  {camera.vp_zv} ft')
print('')
print('11.2.1:')
np.set_printoptions(precision=10, suppress=True)
print(f'x_f_vp:  {camera.vp_xf}')
np.set_printoptions(precision=11, suppress=True)
print(f'y_f_vp:  {camera.vp_yf}')
print(f'z_f_vp:  {camera.vp_zf}')
print('')
print('11.3.1:')
print(f'P0:      {camera.P0}')
print(f'n_vp:    {camera.nvp}')
print(f'P0 - PC: {(camera.P0 - camera_location.flatten())}')
print('lca: ')
# print(ground.lca)
# print('gamma: ')
# print(ground.gamma)
# np.set_printoptions(precision=12, suppress=True)
# print('x_vp: ')
# print(ground.x_vp)
# print('y_vp: ')
# print(ground.y_vp)

# EXAMPLE SIMULATION
location = camera.location 
quat = camera.quat 

while(location[0] < 0):
    time_begin = time.time()
    camera.set_state(location=location, quat=quat)
    ground.plot(camera=camera)
    fig.canvas.draw()
    fig.canvas.flush_events()
    location[0] += 0.01
    time_end = time.time()
    fps = 1/(time_end - time_begin)
    # print('graphics rate [hz]:', fps)

