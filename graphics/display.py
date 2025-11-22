import json
import numpy as np
import matplotlib.pyplot as plt
import time 

from connection_m import connection
from quaternion import quat_mult, quat_norm, quat_to_euler, eul_to_quat

# OPEN THE JSON FILE
with open('json/graphics.json', 'r') as file:
    graphics = json.load(file)

# CREATE CONNECTION
states_connection = connection(graphics["connections"]["receive_states"])

altitude    = graphics["ground"]["altitude[ft]"]
grid_number = graphics["ground"]["grid_number"]
grid_scale  = graphics["ground"]["grid_scale[ft]"]
color       = graphics["ground"]["color"]

class Camera:
    def __init__(
            self, 
            camera_dict: dict
            ):
        
        # FIND THE FOLLOW DISTANCE
        self.type = camera_dict['type']
        if (self.type == 'follow'):
            self.follow_distance = camera_dict["follow_distance[ft]"]
        else:
            self.follow_distance = 0.0
        # FIND THE VIEW PLANE PROPERTIES
        self.d_vp     = camera_dict['view_plane']['distance[ft]']
        self.theta_vp = np.radians(camera_dict['view_plane']['angle[deg]'])
        self.RA_vp    = camera_dict['view_plane']['aspect_ratio']

        self.wvp = 2.0 * self.d_vp * np.tan(self.theta_vp)
        self.hvp = self.wvp / self.RA_vp

        # FIND THE CAMERA LOCATION AND ORIENTATION
        self.orientation = camera_dict['orientation[deg]']
        self.location = camera_dict['location[ft]']

        # INITIALIZE VIEW PLANE BODY FIXED COORDINATES
        self.xc_vp    =  np.full(4, self.d_vp)

        self.yc_vp    =  np.full(4, self.wvp) * 0.5
        self.yc_vp[0] = -1 * self.yc_vp[0]
        self.yc_vp[1] = -1 * self.yc_vp[1]

        self.zc_vp    =  np.full(4, self.hvp) * 0.5 
        self.zc_vp[0] = -1 * self.zc_vp[0]
        self.zc_vp[3] = -1 * self.zc_vp[3]

        self.camera_fixed = np.array([self.xc_vp, self.yc_vp, self.zc_vp])

        # print(self.camera_fixed)

    def set_state(
            self, 
            location:np.ndarray, 
            quat:np.ndarray
            ):
        
        # SET THE LOCATION AND ORIENTATION 
        self.location = location 
        self.quat = quat 

        # INITIALIZE QUATERNION VALUES
        e0 = self.quat[0]
        ex = self.quat[1]
        ey = self.quat[2]
        ez = self.quat[3]

        self.camera_fixed = np.array([self.xc_vp, self.yc_vp, self.zc_vp])

        # BUILD ROTATION MATRIX AND INVERSE
        rotation_matrix = np.array([[ex**2 + e0**2 - ey**2 - ez**2, 2*(ex*ey + ez*e0), 2*(ex*ez - ey*e0)],
                                    [2*(ex*ey - ez*e0), ey**2 + e0**2 - ex**2 - ez**2, 2*(ey*ez + ex*e0)],
                                    [2*(ex*ez + ey*e0), 2*(ey*ez - ex*e0), ez**2 + e0**2 - ex**2 - ey**2]])
        
        rotation_matrix_inv = np.linalg.inv(rotation_matrix)

        # UPDATE IF FOLLOWING AIRCRAFT
        if(self.type == 'follow'):
            self.location = self.location + np.matmul(rotation_matrix_inv, self.follow_distance)

        # CALCULATE VIEW IN EARTH FIXED PLANE
        self.location = np.array(self.location)
        earth_fixed = np.matmul(rotation_matrix_inv, self.camera_fixed) + self.location[:,np.newaxis]

        self.xf_vp = earth_fixed[0,:]
        self.yf_vp = earth_fixed[1,:]
        self.zf_vp = earth_fixed[2,:]
        print(self.xf_vp)
        print(self.yf_vp)
        print(self.zf_vp)

        # BUILD CENTER POINT OF PLANE
        self.P0 = np.array([[np.mean(earth_fixed[0])], [np.mean(earth_fixed[1])], [np.mean(earth_fixed[2])]]).flatten()

        # CALCULATE THE NORMAL VECTOR 
        P1 = np.array([earth_fixed[0,0], earth_fixed[1,0], earth_fixed[2,0]])
        P2 = np.array([earth_fixed[0,1], earth_fixed[1,1], earth_fixed[2,1]])
        self.nvp = np.cross(P1 - self.P0, P2 - self.P0)
            
class LinesObject:
    def __init__(
            self, 
            lines_dict:dict, 
            ax
            ):
        # SET UP AXES FOR DRAWING
        color = lines_dict["color"]            
        self.ax, = ax.plot([], [], marker='', ls='-', color=color)
        
        # BUILD GROUND PLANE IF SPECIFIED 
        if(lines_dict["type"] == "ground_plane"): 
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

        # BUILD FROM VTK FILE
        if(lines_dict["type"] == "vtk"):
            i 
            
        # ALLOCATE SPACE FOR 3D POINTS IN EARTH FIXED COORDINATES
        self.points3D = np.zeros((self.npoints,3))

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
camera_dict = graphics['camera']
camera = Camera(
            camera_dict=camera_dict
            )
quat = eul_to_quat(camera.orientation)
camera.set_state(location=camera.location, quat=quat)

# PLOT THE FIGURE
fig = plt.figure(figsize=(camera.RA_vp*5.0,5.0))
ax = fig.add_subplot(111)
plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
plt.axis('off')
ax.axes.set_xlim(camera.yc_vp[0], camera.yc_vp[2])
ax.axes.set_ylim(-camera.zc_vp[1], -camera.zc_vp[0])
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.set_xticks([])
ax.set_yticks([])
ax.axes.set_aspect('equal')
fig.canvas.draw()

# BUILD THE GROUND
ground_dict = graphics['ground']
ground = LinesObject(ground_dict, ax)
ground.plot(camera)
plt.show()
# plt.show(block=False)

# # BUILD THE LINES
# ground = Lines(
#             grid_number=grid_number,
#             grid_scale=grid_scale,
#             altitude=altitude,
#             ax=ax
#             )

# location = camera.location 
# quat = camera.quat 

# frame = 0
# fps = 0.0
# states = [0]*14 

# while(frame< 1000):
#     time_begin = time.time()
#     states = states_connection.recv()
#     location = states[7:10]
#     location = np.array(states[7:10])
#     location[2] = -location[2]
#     quat = np.array(states[10:14])
#     # print(location)
#     camera.set_state(location=location, quat=quat)
#     ground.plot(camera=camera)
#     fig.canvas.draw()
#     fig.canvas.flush_events()
#     time_end = time.time()
#     fps = 1/(time_end - time_begin)
#     # print('graphics rate [hz]:', fps)
