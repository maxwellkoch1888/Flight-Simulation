import json
import numpy as np
import matplotlib.pyplot as plt
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
    def __init__(self, grid_number:int, grid_scale:float):
        # BUILD POINTS TO BE PLOTTED 
        xf = []
        yf = []
        zf = []
        self.max_distance = grid_number * grid_scale
        i = -self.max_distance
        while i <= self.max_distance:
            j = -self.max_distance
            while j <= self.max_distance:
                xf.append(i)
                yf.append(j)
                zf.append(altitude)
                j += grid_scale
            i += grid_scale

        self.xf = np.array(xf)
        self.yf = np.array(yf)
        self.zf = np.array(zf)

    def intersecting_points(self, nvp, P0, quat, location):
        # CALCULATE THE LINE FORM THE CAMERA TO EACH POINT OF INTEREST
        lca = [] 
        gamma = []
        point_vp = []
        xvp = []
        yvp = []
        points = np.size(self.xf)
        i = 0

        while i < points:
            point = np.array([[self.xf[i]], [self.yf[i]], [self.zf[i]]])
            distance = point.flatten() - location
            lca.append(distance.flatten())
            # print('point:', point)
            # print('distance:', distance)
            # print('lca:', lca )
            # print('P0:', P0)
            # print('camera_location:', location)

            # CALCULATE GAMMA
            gam = np.matmul((P0-location), nvp) / np.matmul(distance, nvp)
            gamma.append(gam.flatten())

            e0 = quat[0]
            ex = quat[1]
            ey = quat[2]
            ez = quat[3]            
            
            rotation_matrix = np.array([[ex**2 + e0**2 - ey**2 - ez**2, 2*(ex*ey + ez*e0), 2*(ex*ez - ey*e0)],
                                        [2*(ex*ey - ez*e0), ey**2 + e0**2 - ex**2 - ez**2, 2*(ey*ez + ex*e0)],
                                        [2*(ex*ez + ey*e0), 2*(ey*ez - ex*e0), ez**2 + e0**2 - ex**2 - ey**2]])
            
            # CALCULATE xvp, yvp, zvp
            point_vp.append(np.matmul(rotation_matrix, (gam*distance)))
            i += 1

        self.lca = np.array(lca)
        self.gamma = np.array(gamma)
        self.point_vp = np.array(point_vp)
        self.x_vp = self.point_vp[:,1]
        self.y_vp = -self.point_vp[:,2]

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

# BUILD THE LINES
lines = Lines(
            grid_number=grid_number,
            grid_scale=grid_scale
            )
lines.intersecting_points(
    nvp=camera.nvp, 
    P0=camera.P0, 
    quat=quat_init,
    location=camera_location.flatten()
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
print(lines.lca)
print('gamma: ')
print(lines.gamma)
np.set_printoptions(precision=12, suppress=True)
print('x_vp: ')
print(lines.x_vp)
print('y_vp: ')
print(lines.y_vp)

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
# plt.plot(lines.x_vp, lines.y_vp, 'o', color=color) 
# plt.show()
plt.show(block=False)


# # PLOT THE FIGURE
# fig = plt.figure(figsize=(camera.vp_aspect_ratio*5.0,5.0))
# ax = fig.add_subplot(111)
# plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
# plt.axis('off')
# ax.axes.set_xlim(ycvp[0], ycvp[2])
# ax.axes.set_ylim(-zcvp[1], -zcvp[0])
# ax.axes.xaxis.set_ticklabels([])
# ax.axes.yaxis.set_ticklabels([])
# ax.set_xticks([])
# ax.set_yticks([])
# ax.axes.set_aspect('equal')
# fig.canvas.draw()
# # # PLOT POINTS
# # plt.plot(x_vp, y_vp, 'o', color=color) 

# # PLOT HORIZONTAL LINES
# number_of_lines = max_distance * 2 / grid_scale
# nx = int(number_of_lines + 1) # NUMBER OF X POINTS
# ny = int(number_of_lines + 1) # NUMBER OF Y POINTS

# for line_number in range(ny):
#     start = line_number * nx
#     end = start + nx
#     plt.plot(x_vp[start:end], y_vp[start:end], color=color)

# # PLOT VERTICAL LINES
# for line_number in range(nx):
#     start = line_number 
#     end = start + nx * (ny - 1)
#     step = nx
#     plt.plot(x_vp[start::step], y_vp[start::step], color=color)

# plt.show(block=False)
