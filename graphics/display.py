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

        temp_y_coordinate = self.vp_distance * np.tan(0.5*self.vp_angle)
        temp_z_coordinate = temp_y_coordinate / self.vp_aspect_ratio

        # INITIALIZE CAMERA'S BODY FIXED COORDINATES
        self.vp_xv = np.array([self.vp_distance, self.vp_distance, self.vp_distance, self.vp_distance])
        self.vp_yv = np.array([-temp_y_coordinate, -temp_y_coordinate, -temp_y_coordinate, -temp_y_coordinate])
        self.vp_zv = np.array([-temp_z_coordinate, -temp_z_coordinate, -temp_z_coordinate, -temp_z_coordinate])

        # INITIALIZE CAMERA'S EARTH FIXED COORDINATES 
        self.vp_xf = np.array([0.0, 0.0, 0.0, 0.0])
        self.vp_yf = np.array([0.0, 0.0, 0.0, 0.0])
        self.vp_zf = np.array([0.0, 0.0, 0.0, 0.0])
        
        self.lines2D = np.zeros((2,2))

        self.dx = temp_y_coordinate
        self.dy = temp_z_coordinate

        def set_state(self, location:np.ndarray, quat:np.ndarray):
            # SET THE LOCATION AND ORIENTATION 
            self.location = location 
            self.quat = quat 

            


# CALCULATE wvp AND hvp
wvp = 2.0 * dvp * np.tan(theta)
hvp = wvp / RAvp

# BUILD ARRAYS FOR VIEW PLANE
xcvp = np.array([ dvp,  dvp, dvp,  dvp], dtype=np.float64)
ycvp = np.array([-wvp, -wvp, wvp,  wvp], dtype=np.float64) * 0.5
zcvp = np.array([-hvp,  hvp, hvp, -hvp], dtype=np.float64) * 0.5

# CALCULATE VIEW IN EARTH FIXED PLANE
e0 = camera_orientation_quat[0]
ex = camera_orientation_quat[1]
ey = camera_orientation_quat[2]
ez = camera_orientation_quat[3]
camera_fixed = np.array([xcvp, ycvp, zcvp])

rotation_matrix = np.array([[ex**2 + e0**2 - ey**2 - ez**2, 2*(ex*ey + ez*e0), 2*(ex*ez - ey*e0)],
                            [2*(ex*ey - ez*e0), ey**2 + e0**2 - ex**2 - ez**2, 2*(ey*ez + ex*e0)],
                            [2*(ex*ez + ey*e0), 2*(ey*ez - ex*e0), ez**2 + e0**2 - ex**2 - ey**2]])

earth_fixed = np.matmul(np.linalg.inv(rotation_matrix), camera_fixed) + camera_location

# BUILD POINTS TO BE PLOTTED 
xf = []
yf = []
zf = []
max_distance = grid_number * grid_scale
i = -max_distance
while i <= max_distance:
    j = -max_distance
    while j <= max_distance:
        xf.append(i)
        yf.append(j)
        zf.append(altitude)
        j += grid_scale
    i += grid_scale

xf = np.array(xf)
yf = np.array(yf)
zf = np.array(zf)

# BUILD CENTER POINT OF PLANE
P0 = np.array([[np.mean(earth_fixed[0])], [np.mean(earth_fixed[1])], [np.mean(earth_fixed[2])]])

# CALCULATE THE NORMAL VECTOR 
P1 = np.array([earth_fixed[0,0], earth_fixed[1,0], earth_fixed[2,0]])
P2 = np.array([earth_fixed[0,1], earth_fixed[1,1], earth_fixed[2,1]])
P0 = P0.flatten()
camera_location = camera_location.flatten()

nvp = np.cross(P1 - P0, P2 - P0)

# CALCULATE THE LINE FORM THE CAMERA TO EACH POINT OF INTEREST
lca = [] 
gamma = []
point_vp = []
xvp = []
yvp = []
points = np.size(xf)
i = 0

while i < points:
    point = np.array([[xf[i]], [yf[i]], [zf[i]]])
    distance = point.flatten() - camera_location
    lca.append(distance.flatten())

    # CALCULATE GAMMA
    gam = np.matmul((P0-camera_location), nvp) / np.matmul(distance, nvp)
    gamma.append(gam.flatten())

    # CALCULATE xvp, yvp, zvp
    point_vp.append(np.matmul(rotation_matrix, (gam*distance)))
    i += 1
lca = np.array(lca)
gamma = np.array(gamma)
point_vp = np.array(point_vp)
x_vp = point_vp[:,1]
y_vp = -point_vp[:,2]



# PRINT RESULTS
print('------------ Results ------------')
print(f'w_vp:    {wvp:.11f} ft')
print(f'h_vp:    {hvp:.11f} ft')
print('')
np.set_printoptions(precision=2, suppress=True)
print(f'x_c_vp:  {xcvp} ft')
np.set_printoptions(precision=11, suppress=True)
print(f'y_c_vp:  {ycvp} ft')
print(f'z_c_vp:  {zcvp} ft')
print('')
np.set_printoptions(precision=10, suppress=True)
print(f'x_f_vp:  {earth_fixed[0]}')
np.set_printoptions(precision=11, suppress=True)
print(f'y_f_vp:  {earth_fixed[1]}')
print(f'z_f_vp:  {earth_fixed[2]}')
print('')
print(f'P0:      {P0.flatten()}')
print(f'n_vp:    {nvp}')
print(f'P0 - PC: {(P0 - camera_location).flatten()}')
print('lca: ')
print(lca)
print('gamma: ')
print(gamma)
np.set_printoptions(precision=12, suppress=True)
print('x_vp: ')
print(x_vp)
print('y_vp: ')
print(y_vp)

# BUILD THE CAMERA
camera = Camera(
            orientation= camera_orientation,
            location=camera_location, 
            vp_distance=dvp,
            vp_angle=theta,
            vp_aspect_ratio=RAvp
            )

# PLOT THE FIGURE
fig = plt.figure(figsize=(camera.vp_aspect_ratio*5.0,5.0))
ax = fig.add_subplot(111)
plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
plt.axis('off')
ax.axes.set_xlim(ycvp[0], ycvp[2])
ax.axes.set_ylim(-zcvp[1], -zcvp[0])
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.set_xticks([])
ax.set_yticks([])
ax.axes.set_aspect('equal')
fig.canvas.draw()
# # PLOT POINTS
# plt.plot(x_vp, y_vp, 'o', color=color) 

# PLOT HORIZONTAL LINES
number_of_lines = max_distance * 2 / grid_scale
nx = int(number_of_lines + 1) # NUMBER OF X POINTS
ny = int(number_of_lines + 1) # NUMBER OF Y POINTS

for line_number in range(ny):
    start = line_number * nx
    end = start + nx
    plt.plot(x_vp[start:end], y_vp[start:end], color=color)

# PLOT VERTICAL LINES
for line_number in range(nx):
    start = line_number 
    end = start + nx * (ny - 1)
    step = nx
    plt.plot(x_vp[start::step], y_vp[start::step], color=color)

plt.show()
